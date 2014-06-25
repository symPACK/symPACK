all: ETree.o cost.o
	g++ -I./ -o md main.cpp cost.o ETree.o -O0 -g

#	g++ -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o: ETree.cpp
	g++ -I./ -o ETree.o -c ETree.cpp -O0 -g -std=c++11

cost.o: cost.cpp
	g++ -I./ -o cost.o -c cost.cpp -O0 -g -std=c++11

call_cost: ETree.o cost.o
	g++ -o cost call_cost.cpp ETree.o cost.o -O0 -g

clean:
	rm -rf cost md *.o
