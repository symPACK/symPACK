from random import shuffle
import sys

num_genes = int(sys.argv[1])
pop_size = int(sys.argv[2])
i = list(range(1, num_genes + 1))
shuffle(i)
for j in range(0,pop_size):
	print(*i, sep=' ')
	shuffle(i)