CXX = g++ -O3 -std=c++17 
CXXFLAGS = -fopenmp -pthread -march=native -mtune=native -O3 -lm -lz -lpthread -ldl

all: DipGenie

src_dir := src

OBJS = $(src_dir)/main.o $(src_dir)/gfa-io.o $(src_dir)/gfa-base.o \
		$(src_dir)/options.o $(src_dir)/kalloc.o \
		$(src_dir)/sys.o $(src_dir)/approximator.o $(src_dir)/MurmurHash3.o \
		$(src_dir)/misc.o $(src_dir)/solver.o

DipGenie: $(OBJS)
	$(CXX) $^ -o $@ $(CXXFLAGS)

$(src_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(src_dir)/*.o DipGenie
