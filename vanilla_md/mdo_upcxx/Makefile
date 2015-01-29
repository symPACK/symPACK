EXE = mdo_upcxx
#EXE = mdo_explicit
#EXE = mdo_reach
#EXE = mdo_explicit_queue

.PHONY: all clean 

all: clean mdo_reach mdo_explicit mdo_explicit_queue
	make $(EXE)

clean:
	-rm -f *.o $(EXE)

#UPCXX_INSTALL_DIR ?= $(HOME)/upcxx_install
#UPCXX_INSTALL_DIR ?= $(HOME)/install/upcxx_clang
#UPCXX_INSTALL_DIR ?= $(PREFIX)/upcxx_install
UPCXX_INSTALL_DIR ?= $(HOME)/software/release/upcxx_install
include $(UPCXX_INSTALL_DIR)/include/upcxx.mak

VPATH=$(UPCXX_INSTALL_DIR)/include

MYCXXFLAGS =-O3 -fopenmp -DPROFILE
#MYCXXFLAGS =-O0 -g -fopenmp -DPROFILE
CXX=mpicxx 


util.o: util.cpp
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

ETree.o: ETree.cpp
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

timer.o: timer.cpp
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

$(EXE).o: $(EXE).cpp util.cpp ETree.cpp timer.cpp 
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

$(EXE): $(EXE).o util.o ETree.o timer.o
	$(CXX) $(LDFLAGS) $(UPCXX_LDFLAGS) $(MYCXXFLAGS) -o $@ $^ $(UPCXX_LDLIBS) 

mdo_explicit.o: mdo_explicit.cpp util.cpp ETree.cpp timer.cpp 
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

mdo_explicit: mdo_explicit.o util.o ETree.o timer.o
	$(CXX) $(LDFLAGS) $(UPCXX_LDFLAGS) $(MYCXXFLAGS) -o $@ $^ $(UPCXX_LDLIBS) 

mdo_explicit_queue.o: mdo_explicit_queue.cpp util.cpp ETree.cpp timer.cpp 
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

mdo_explicit_queue: mdo_explicit_queue.o util.o ETree.o timer.o
	$(CXX) $(LDFLAGS) $(UPCXX_LDFLAGS) $(MYCXXFLAGS) -o $@ $^ $(UPCXX_LDLIBS) 

mdo_reach.o: mdo_reach.cpp util.cpp ETree.cpp timer.cpp 
	$(CXX) $(UPCXX_CXXFLAGS) -I. $(MYCXXFLAGS) -o $@ -c $<

mdo_reach: mdo_reach.o util.o ETree.o timer.o
	$(CXX) $(LDFLAGS) $(UPCXX_LDFLAGS) $(MYCXXFLAGS) -o $@ $^ $(UPCXX_LDLIBS) 
