#DEBUG_LIB = ${HOME}/software/duma_2_5_15/libduma.a -lpthread

CXX=CC
CXXFLAGS =-I./  -O3 -std=c++11
#CXXFLAGS =-I./  -g -O0 -std=c++11

all: check_md md vmd call_cost maxd genetic

vmd: ETree.o util.o main.cpp 
	${CXX} ${CXXFLAGS} -o vmd main.cpp util.o ETree.o 

genetic: geneticOrdering.c util.o ETree.o 
	${CXX} ${CXXFLAGS} -fopenmp -o geneticOrdering geneticOrdering.c util.o ETree.o

md: minimum_degree.cpp util.o ETree.o
	${CXX} ${CXXFLAGS} -o md minimum_degree.cpp util.o ETree.o

maxd: maximum_degree.cpp util.o ETree.o
	${CXX} ${CXXFLAGS} -o maxd maximum_degree.cpp util.o ETree.o
#	${CXX} -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o: ETree.cpp
	${CXX} ${CXXFLAGS} -o ETree.o -c ETree.cpp 

util.o: util.cpp ETree.o
	${CXX} ${CXXFLAGS} -o util.o -c util.cpp 

call_cost: ETree.o util.o call_cost.cpp
	${CXX} ${CXXFLAGS} -o call_cost call_cost.cpp ETree.o util.o ${DEBUG_LIB}


fastshift: ETree.o util.o fastshift.cpp
	${CXX} ${CXXFLAGS} -o fastshift fastshift.cpp ETree.o util.o 


check_md: ETree.o util.o check_md.cpp
	${CXX} ${CXXFLAGS} -o check_md check_md.cpp ETree.o util.o 

clean:
	rm -rf core call_cost check_md md mymd *.o geneticOrdering

