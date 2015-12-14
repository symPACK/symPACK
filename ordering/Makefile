FC=gfortran

all: box grid

box: box_main.f boxnd.f
	$(FC) -o box box_main.f boxnd.f

grid: grid_main.f gridnd.f
	$(FC) -o grid grid_main.f gridnd.f
