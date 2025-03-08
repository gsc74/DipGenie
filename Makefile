CXX = g++ -std=c++17
CXXFLAGS = -fopenmp -pthread -march=native -mtune=native -O3
LDFLAGS = -lgurobi_c++ -lgurobi120 -lm -lz -lpthread -ldl
GUROBI_HOME = /home/ghanshyam/opt/gurobi1201/linux64
INCLUDES = -Iextra/include
LIBS = -Lextra/lib

src_dir := src

OBJS = $(src_dir)/main.o $(src_dir)/gfa-io.o $(src_dir)/gfa-base.o \
       $(src_dir)/options.o $(src_dir)/kalloc.o \
       $(src_dir)/misc.o $(src_dir)/sys.o $(src_dir)/ILP_index.o \
       $(src_dir)/MurmurHash3.o

all: PHI2

PHI2: $(OBJS)
	$(CXX) $^ -o $@ $(INCLUDES) $(LIBS) $(CXXFLAGS) $(LDFLAGS)

$(src_dir)/%.o: $(src_dir)/%.cpp
	$(CXX) -c $< -o $@ $(INCLUDES) $(LIBS) $(CXXFLAGS)

clean:
	rm -f $(src_dir)/*.o PHI2