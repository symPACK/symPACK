all: ETree.o
	g++ -o md main.cpp ETree.o -O0 -g
	g++ -o md_set main_set.cpp -O0 -g -std=c++11

ETree.o:
	g++ -I./ -o ETree.o -c ETree.cpp -O0 -g -std=c++11
