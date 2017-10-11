exe = subgraph 

cc = mpiCC#"$(shell which g++)" 
flags = -I. -fopenmp -march=athlon64 -O3 -ggdb
#flags = -I. -fopenmp
flags += -g -gdwarf-2
flags += -std=c++11

#ifeq ($(debug), 1)
#	flags+= -O0 -g #-DDEBUG 
#else
#	flags += -O3
#endif

objs = $(patsubst %.cpp,%.o,$(wildcard ../../lib/*.cpp)) \
			$(patsubst %.cpp,%.o,$(wildcard *.cpp))

deps = $(wildcard ./*.hpp) \
				$(wildcard *.h) \
				Makefile

%.o:%.cpp $(deps)
	$(cc) -c $< -o $@ $(flags)

$(exe):$(objs)
	$(cc) $(objs) -o $(exe) $(flags)

clean:
	rm -rf $(exe) $(objs) 
