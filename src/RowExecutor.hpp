#pragma once
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <semaphore.h>
#include <stdexcept>
#include <thread>
#include <vector>

namespace hic {
// Persistent, sleeping workers for independent destination rows. The caller
// participates; one completion counter protects each row batch's lifetime.
// No busy-wait barrier and no order-dependent reduction are used.
class RowExecutor {
    sem_t ready_;
    std::vector<std::thread> workers_;
    std::atomic<bool> stop_{false},failed_{false};
    std::atomic<size_t> next_{0},pending_{0},epoch_{0};
    size_t count_=0;
    std::function<void(size_t)> task_;
    std::exception_ptr failure_;
    std::mutex finish_mutex_,failure_mutex_;
    std::condition_variable finished_;

    void consume() noexcept {
        try {
            while(!failed_.load(std::memory_order_relaxed)) {
                const size_t i=next_.fetch_add(1,std::memory_order_relaxed);
                if(i>=count_)break;
                task_(i);
            }
        } catch(...) {
            std::lock_guard<std::mutex> lock(failure_mutex_);
            if(!failure_)failure_=std::current_exception();
            failed_.store(true,std::memory_order_relaxed);
        }
        if(pending_.fetch_sub(1,std::memory_order_acq_rel)==1) {
            std::lock_guard<std::mutex> lock(finish_mutex_);
            finished_.notify_one();
        }
    }
    void worker() noexcept {
        for(;;) {
            int code;
            do {code=sem_wait(&ready_);}while(code<0&&errno==EINTR);
            if(code<0)std::terminate();
            if(stop_.load(std::memory_order_acquire))return;
            (void)epoch_.load(std::memory_order_acquire);
            consume();
        }
    }
    void join() noexcept {
        stop_.store(true,std::memory_order_release);
        for(size_t i=0;i<workers_.size();++i)if(sem_post(&ready_)<0)std::terminate();
        for(auto&t:workers_)t.join();
    }
public:
    explicit RowExecutor(int threads) {
        if(threads<1)throw std::invalid_argument("positive row-executor thread count required");
        if(sem_init(&ready_,0,0)<0)throw std::runtime_error("Cannot initialize row semaphore");
        try {for(int i=1;i<threads;++i)workers_.emplace_back([this]{worker();});}
        catch(...){join();sem_destroy(&ready_);throw;}
    }
    ~RowExecutor(){join();sem_destroy(&ready_);}
    RowExecutor(const RowExecutor&)=delete;
    RowExecutor&operator=(const RowExecutor&)=delete;

    template<class F> void for_each(size_t count,F&&function) {
        if(!count)return;
        if(workers_.empty()||count==1){for(size_t i=0;i<count;++i)function(i);return;}
        task_=std::forward<F>(function);count_=count;failure_=nullptr;
        next_.store(0,std::memory_order_relaxed);failed_.store(false,std::memory_order_relaxed);
        const size_t helpers=std::min(workers_.size(),count-1);
        pending_.store(helpers+1,std::memory_order_relaxed);
        epoch_.fetch_add(1,std::memory_order_release);
        for(size_t i=0;i<helpers;++i)if(sem_post(&ready_)<0)std::terminate();
        consume();
        {
            std::unique_lock<std::mutex> lock(finish_mutex_);
            finished_.wait(lock,[this]{return pending_.load(std::memory_order_acquire)==0;});
        }
        task_={};
        if(failure_)std::rethrow_exception(failure_);
    }
};
}
