#!/usr/bin/python
from __future__ import print_function
import sys
import re
import fileinput
import numpy
import collections
import csv
import math
import numbers

def lookahead(iterable):
    it = iter(iterable)
    last = next(it) # next(it) in Python 3
    for val in it:
        yield last, False
        last = val
    yield last, True

def nplookahead(iterable):
    it = np.nditer(iterable)
    last = it;
    last.iternext() # next(it) in Python 3
    
    for val in it:
        yield last, False
        last = val
    yield last, True


def flatten(d, level=-1, parent_key='', sep=' '):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            if(level == -1):
              items.extend(flatten(v,level, new_key, sep=sep).items())
            elif(level>0):
              items.extend(flatten(v,level-1, new_key, sep=sep).items())
            else:
              items.append((new_key, v))
        else:
            items.append((new_key, v))
    return dict(items)


matrixre = re.compile("[ -]+/data/matrices/([a-zA-Z_0-9-]+)/")
orderre = re.compile("[ ]*[-]+ ([a-zA-Z]+) [-]+")

sblockre = re.compile("[ ]*totalSnodeBlocks:[ ]*([0-9.]+)")
blockre = re.compile("[ ]*totalBlocks:[ ]*([0-9.]+)")
timere = re.compile("[ ]*Refinement done in ([0-9.]+)")

matrix=""
order=""
matrices={}

#type%2==1
# 0 = TSP, 1 = Barry 2 1, 2 = Barry 2 0, 3 = Barry 1 1, 4 = Barry 1 0
type = -1
ttype = -1
for line in fileinput.input():
    #print line
    matrixm = matrixre.search(line)
    orderm = orderre.search(line)

    if(matrixm):
      matrix = matrixm.group(1)

    if(orderm):
      order = orderm.group(1)

      ttype = -1
      type=-1

    
    sblockm = sblockre.search(line)
    blockm = blockre.search(line)
    timem = timere.search(line)

    if(timem):
      ttype = ttype+1
      time = timem.group(1)

      if(ttype == 0):
        mtype = 'TSP'
      elif(ttype == 1):
        mtype = 'Barry_2_1'
      elif(ttype == 2):
        mtype = 'Barry_2_0'
      elif(ttype == 3):
        mtype = 'Barry_1_1'
      elif(ttype == 4):
        mtype = 'Barry_1_0'


      #do the insertion
      if(not order in matrices):
          matrices[order] = {}

      if(not "time" in matrices[order]):
        matrices[order]["time"] = {}
      if(not mtype in matrices[order]["time"]):
          matrices[order]["time"][mtype] = {}
      if(not matrix in matrices[order]["time"][mtype]):
          matrices[order]["time"][mtype][matrix] = []
      matrices[order]["time"][mtype][matrix].append(float(time))





    if(sblockm):
      type = type+1
      mtype = str(type/2)

      if(type/2 == 0):
        mtype = 'TSP'
      elif(type/2 == 1):
        mtype = 'Barry_2_1'
      elif(type/2 == 2):
        mtype = 'Barry_2_0'
      elif(type/2 == 3):
        mtype = 'Barry_1_1'
      elif(type/2 == 4):
        mtype = 'Barry_1_0'


      sblock = sblockm.group(1)

      #do the insertion
      if(not order in matrices):
          matrices[order] = {}

      if(type == 0):
        if(not "No refinement" in matrices[order]):
          matrices[order]["No refinement"] = {}
        if(not mtype in matrices[order]["No refinement"]):
            matrices[order]["No refinement"][mtype] = {}
        if(not matrix in matrices[order]["No refinement"][mtype]):
            matrices[order]["No refinement"][mtype][matrix] = []
        matrices[order]["No refinement"][mtype][matrix].append(int(sblock))
      elif(type % 2 == 1):
        if(not "superBlocks" in matrices[order]):
          matrices[order]["superBlocks"] = {}
        if(not mtype in matrices[order]["superBlocks"]):
            matrices[order]["superBlocks"][mtype] = {}
        if(not matrix in matrices[order]["superBlocks"][mtype]):
            matrices[order]["superBlocks"][mtype][matrix] = []
        matrices[order]["superBlocks"][mtype][matrix].append(int(sblock))


#    if(blockm):
#      block = blockm.group(1)
#      #do the insertion
#      if(not order in matrices):
#          matrices[order] = {}
#
#      if(not "blocks" in matrices[order]):
#        matrices[order]["blocks"] = {}
#
#      if(not type in matrices[order]["blocks"]):
#          matrices[order]["blocks"][type] = {}
#
#      if(not matrix in matrices[order]["blocks"][type]):
#          matrices[order]["blocks"][type][matrix] = []
#
#      matrices[order]["blocks"][type][matrix].append(int(block))
#print(matrices)
flat = flatten(matrices,2)


headers = []

headers = flat.keys()

#it =  iter( flat.keys())
#for key in it:
#  if(len(flat[key].keys())>len(headers)):
#    headers = flat[key].keys()


matrixnames = []
print("%s, " % ('matrix'), end='')
for header,last in lookahead(headers):
  if(not(last)):
    print("%s, " % (header), end='')
  else:
    print("%s" % (header))

  if(len(flat[header].keys())>len(matrixnames)):
     matrixnames = flat[header].keys()


#print(flatten(matrices,1))
#print(flatten(matrices,2))

for matrix in matrixnames:
  print("%s, " % (matrix), end='')
  for header,last in lookahead(flat.keys()):
    val = 'nan'
    if matrix in flat[header]:
      if(len(flat[header][matrix])>0):
        if(isinstance(flat[header][matrix][0],numbers.Integral)):
          val = ('%d' % math.ceil(numpy.nanmean(flat[header][matrix])))
        else:
          val = ('%.5f' % numpy.nanmean(flat[header][matrix]))

    if(not(last)):
      print("%s, " % (("%s" % val)), end='')
    else:
      print("%s" % (("%s" % val)))

#  print(header)

#  maxVals=[]
#  for header,last in lookahead(flat[head].keys()):
#    if(len(flat[head][header])>len(maxVals)):
#      maxVals = flat[head][header]
#
#    print(maxVals)
#  for i in range(0,len(maxVals)):  
#    for type, last in lookahead(headers):
#      val = float('nan')
#      if(type in flat[head].keys()):
#        if(i<len(flat[head][type])):
#          val = valuesDict[proc][type][i]
#  #    #print(val)
#      if(not(last)):
#        print("%s, " % (("%.5f" % val).rjust(10,' ')), end='')
#      else:
#        print("%s" % (("%.5f" % val).rjust(10,' ')))


