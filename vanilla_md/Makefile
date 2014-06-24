all: ETree.o cost.o
	g++ -I./ -o md main.cpp cost.o ETree.o -O0 -g

#	g++ -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o:
	g++ -I./ -o ETree.o -c ETree.cpp -O0 -g -std=c++11

cost.o:
	g++ -I./ -o cost.o -c cost.cpp -O0 -g -std=c++11
