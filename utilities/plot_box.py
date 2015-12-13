# !/bin/python

import fileinput
import math
import sys

import numpy as np
from mpl_toolkits import mplot3d
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

k = int(sys.argv[1])

# i is the point we are interested into
i = int(sys.argv[2])-1
points = [int(x)-1 for x in sys.argv[3:]]

# test of a linear to 3d coordinates corresponding to what we'd like
def to3dcoords(i):
  xDirection = (i) % k;
  yDirection = ((i) / k) % k;
  zDirection = (i) / (k**2); 
  #print "%d=(%d,%d,%d)" % (i,xDirection,yDirection,zDirection)
  return [xDirection,yDirection,zDirection]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = [item for sublist in [range(0,k)*k**2] for item in sublist]
ys = [item for sublist in [ [x]*k for x in range(0,k) ]*k for item in sublist]
zs = [item for sublist in [ [x]*k**2 for x in range(0,k) ] for item in sublist]
print xs
print ys
print zs
#xs = range(0,k**3)
#ys = range(0,k**3)
#ys = range(0,k**3)
cs=['r']*k**3
cs[i] = 'g'
for idx in points:
  cs[idx] = 'b'

ax.scatter(xs, ys, zs, c=cs)#, cmap=plt.get_cmap('nipy_spectral'))


#labels = ['{0}'.format(i) for i in range(0,k**3)]


#lines = iter([((1, 1, 1), (2, 2, 2))])
#l = mplot3d.art3d.Line3DCollection(lines)
#ax.add_collection(l)
segments = []
for p in points:
  segments.append([ to3dcoords(i), to3dcoords(p) ])
print segments
#segments [list(zip(x, y)) for y in ys]

l = mplot3d.art3d.Line3DCollection(segments)
#l.set_edgecolor('k')
ax.add_collection3d(l)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
