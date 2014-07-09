CXXFLAGS =-I./ -O3 -std=c++11
#CXXFLAGS =-I./ -O0 -g -std=c++11 


all: check_md md vmd call_cost maxd

vmd: ETree.o util.o main.cpp 
	g++ ${CXXFLAGS} -o vmd main.cpp util.o ETree.o 

genetic: geneticOrdering.c util.o ETree.o 
	g++ ${CXXFLAGS} -o geneticOrdering geneticOrdering.c util.o Etree.o

md: minimum_degree.cpp util.o ETree.o
	g++ ${CXXFLAGS} -o md minimum_degree.cpp util.o ETree.o

maxd: maximum_degree.cpp util.o ETree.o
	g++ ${CXXFLAGS} -o maxd maximum_degree.cpp util.o ETree.o
#	g++ -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o: ETree.cpp
	g++ ${CXXFLAGS} -o ETree.o -c ETree.cpp 

util.o: util.cpp ETree.o
	g++ ${CXXFLAGS} -o util.o -c util.cpp 

call_cost: ETree.o util.o call_cost.cpp
	g++ ${CXXFLAGS} -o call_cost call_cost.cpp ETree.o util.o 


fastshift: ETree.o util.o fastshift.cpp
	g++ ${CXXFLAGS} -o fastshift fastshift.cpp ETree.o util.o 


check_md: ETree.o util.o check_md.cpp
	g++ ${CXXFLAGS} -o check_md check_md.cpp ETree.o util.o 

clean:
	rm -rf call_cost check_md md mymd *.o geneticOrdering

