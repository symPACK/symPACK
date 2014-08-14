#import numpy as np
import matplotlib
matplotlib.use('GTK3Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys

if len(sys.argv) != 2:
    print("Only takes one imput file")
inputfile = str(sys.argv[1])
openedFile = open(inputfile, 'r')

largestX = 0
largestY = 0
dataList = []

nextLine = openedFile.readline()
while nextLine != "":
    splitLine = nextLine.split()
    dataList.append([float(splitLine[1]), float(splitLine[2]), splitLine[3], float(splitLine[4])])
    if float(splitLine[4]) > largestY:
        largestY = float(splitLine[4])
    if float(splitLine[2]) > largestX:
        largestX = float(splitLine[2])
    nextLine = openedFile.readline()

openedFile.close()
#print dataList

#plt.ion()

def printLabel(event):
    print event.artist.get_gid()

fig = plt.figure() 
axes = plt.gca()
axes.set_xlim(0,largestX + .5)
axes.set_ylim(-.5,largestY + .5)

for el in dataList:
    p = mpatches.Rectangle((el[0], el[3]-.25), el[1]-el[0], .5, picker=True, gid=el[2])
    if el[2][0] == 'U':
        p.set_color('r')
    axes.add_patch(p)

fig.canvas.mpl_connect('pick_event', printLabel)
plt.show()

