CXX = g++ -O3 -std=c++17 -DUSE_BAYESOPT
CXXFLAGS = -fopenmp -pthread -march=native -mtune=native -O3 -lgurobi_c++ -lgurobi120 -lm -lz -lpthread -ldl -lbayesopt -lnlopt
GUROBI_HOME = /home/ghanshyam/opt/gurobi1202/linux64
INCLUDES = -I./extra/include/ -I/data/nfs_home/ghanshy1/RECOMB/bayesopt/include
LIBS = -L./extra/lib/ -L/data/nfs_home/ghanshy1/RECOMB/bayesopt/lib

all: PHI2

src_dir := src

OBJS = $(src_dir)/main.o $(src_dir)/gfa-io.o $(src_dir)/gfa-base.o \
		$(src_dir)/options.o $(src_dir)/kalloc.o \
		$(src_dir)/sys.o $(src_dir)/ILP_index.o \
		$(src_dir)/approximator.o $(src_dir)/MurmurHash3.o \
		$(src_dir)/misc.o

PHI2: $(OBJS)
	$(CXX) $^ -o $@ $(INCLUDES) $(LIBS) $(CXXFLAGS)

$(src_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) -c $< -o $@ $(INCLUDES) $(LIBS) $(CXXFLAGS)

clean:
	rm -f $(src_dir)/*.o PHI2
