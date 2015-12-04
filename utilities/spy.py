# !/bin/python
import matplotlib.pylab as plt
import scipy.sparse as sps
import scipy.io as sio
import sys

A = sio.mmread(sys.argv[1])

#A = sps.rand(10000,10000, density=0.00001)
M = sps.csr_matrix(A)
plt.spy(M,markersize=1)
plt.show()
