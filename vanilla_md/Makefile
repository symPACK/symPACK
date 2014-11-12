#DEBUG_LIB = ${HOME}/software/duma_2_5_15/libduma.a -lpthread

FC=ftn
CXX=CC

FFLAGS=

CXXFLAGS =-I./ -DAdd_ -fast -std=c++11 -openmp -DPROFILE
#CXXFLAGS =-I./  -DAdd_ -g -O0 -std=c++11 -openmp 
#-DPROFILE -openmp

all: check_md md vmd call_cost maxd genetic

SRCS_F= genmmd.f  mmdelm.f  mmdint.f  mmdnum.f  mmdupd.f  ordmmd.f 

FORTRANOBJS= ${SRCS_F:.f=.o}


mmd: ${FORTRANOBJS}
	${CXX} ${CXXFLAGS} -o mmd mmd.cpp util.o ETree.o timer.o ${FORTRANOBJS}


%.o: %.f
	${FC} -c ${FFLAGS} $<

vmd: ETree.o util.o main.cpp 
	${CXX} ${CXXFLAGS} -o vmd main.cpp util.o ETree.o 

genetic: geneticOrdering.c util.o ETree.o 
	${CXX} ${CXXFLAGS} -fopenmp -o geneticOrdering geneticOrdering.c util.o ETree.o

md: minimum_degree.cpp util.o ETree.o timer.o
	${CXX} ${CXXFLAGS} -o md minimum_degree.cpp util.o ETree.o timer.o

md_list: minimum_degree_list.cpp util.o ETree.o timer.o
	${CXX} ${CXXFLAGS} -o md_list minimum_degree_list.cpp util.o ETree.o timer.o

maxd: maximum_degree.cpp util.o ETree.o
	${CXX} ${CXXFLAGS} -o maxd maximum_degree.cpp util.o ETree.o
#	${CXX} -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o: ETree.cpp
	${CXX} ${CXXFLAGS} -o ETree.o -c ETree.cpp 

util.o: util.cpp ETree.o
	${CXX} ${CXXFLAGS} -o util.o -c util.cpp 

timer.o: timer.cpp 
	${CXX} ${CXXFLAGS} -o timer.o -c timer.cpp 

call_cost: ETree.o util.o call_cost.cpp
	${CXX} ${CXXFLAGS} -o call_cost call_cost.cpp ETree.o util.o ${DEBUG_LIB}


fastshift: ETree.o util.o fastshift.cpp
	${CXX} ${CXXFLAGS} -o fastshift fastshift.cpp ETree.o util.o 


check_md: ETree.o util.o check_md.cpp
	${CXX} ${CXXFLAGS} -o check_md check_md.cpp ETree.o util.o 

clean:
	rm -rf core call_cost check_md md mymd *.o geneticOrdering

