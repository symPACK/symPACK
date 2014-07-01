#CXXFLAGS =-I./ -O0 -g -std=c++11 
CXXFLAGS =-I./ -O3 -std=c++11 


all: ETree.o cost.o main.cpp md
	g++ ${CXXFLAGS} -o mymd main.cpp cost.o ETree.o 

md: minimum_degree.cpp
	g++ ${CXXFLAGS} -o md minimum_degree.cpp 

maxd: maximum_degree.cpp
	g++ ${CXXFLAGS} -o maxd maximum_degree.cpp 
#	g++ -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o: ETree.cpp
	g++ ${CXXFLAGS} -o ETree.o -c ETree.cpp 

cost.o: cost.cpp
	g++ ${CXXFLAGS} -o cost.o -c cost.cpp 

call_cost: ETree.o cost.o
	g++ ${CXXFLAGS} -o call_cost call_cost.cpp ETree.o cost.o 

check_md: ETree.o cost.o check_md.cpp
	g++ ${CXXFLAGS} -o check_md check_md.cpp ETree.o cost.o 
clean:
	rm -rf call_cost check_md md mymd *.o
