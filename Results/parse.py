#!/usr/bin/python
from __future__ import print_function
import sys
import re
import fileinput
import numpy

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





pullregex = re.compile("run_sparse_pull")
pushregex = re.compile("run_sparse_push")
mumpsregex = re.compile("run_mumps")
pastixregex = re.compile("run_pastix")
elementalregex = re.compile("run_elemental")
superluregex = re.compile("run_superlu")
mapregex = re.compile("-map ([a-zA-Z0-9]+)")

tregex = re.compile("[ ]*Factorization time:[ ]*([0-9.]+)[ ]*[s]?")
tregex2 = re.compile("^Factorization done in: ([0-9.]+)s")
tregex3 = re.compile("[ ]*Time to factorize[ ]*([0-9.]+) s")
tregex4 = re.compile("^Time for factorization is ([0-9.]+) sec")
pregex = re.compile("^[-]+([0-9]+)[-]+$")
pregex2 = re.compile("[ ]*NUMBER OF WORKING PROCESSES[ ]*=[ ]*([0-9]+)")


errorregex = re.compile("ERROR")

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

values={}
valuesDict={}
output=0
type=""
mapping=""
haserror=0
#with open(sys.argv[1], "rtU") as f:
#  for line in f:
for line in fileinput.input():
    #print line
    pull = pullregex.search(line)
    push = pushregex.search(line)
    pastix = pastixregex.search(line)
    superlu = superluregex.search(line)
    elemental = elementalregex.search(line)
    mumps = mumpsregex.search(line)
    map = mapregex.search(line)
    error = errorregex.search(line)

    if(error):
      haserror=1

    if(pull):
      type="(PULL)"
      haserror=0
    elif(push):
      type="(PUSH)"
      haserror=0
    elif(mumps):
      type="\"MUMPS 5.0\""
      haserror=0
    elif(pastix):
      type="\"PASTIX 5.2.2\""
      haserror=0
    elif(superlu):
      type="\"SuperLU_DIST 3.3\""
      haserror=0
    elif(elemental):
      type="\"ELEMENTAL\""
      haserror=0

    if(map):
      mapping = map.group(1)

    m = pregex.match(line)
    m2 = pregex2.match(line)
    if(m):
      proc = int(m.group(1))
#      print "\n%d,"%proc,
    elif(m2):
      proc = int(m2.group(1))
#      print "\n%d,"%proc,
    else:
      m = tregex.match(line)
      m2 = tregex2.match(line)
      m3 = tregex3.match(line)
      m4 = tregex4.match(line)

      time="0"
      match=False
      if(m):
        if(type != "\"ELEMENTAL\""):
          type="\""+mapping+" "+type+"\""
        if(haserror):
          time="nan"
        else:
          time = m.group(1)
        match=True
      elif(m2):
        if(haserror):
          time="nan"
        else:
          time = m2.group(1)
#        print ("%.5f," % float(time)),
        match=True
      elif(m3):
        if(haserror):
          time="nan"
        else:
          time = m3.group(1)
#        print ("%.5f," % float(time)),
        match=True
      elif(m4):
        if(haserror):
          time="nan"
        else:
          time = m4.group(1)
#        print ("%.5f," % float(time)),
        match=True

      if(match):
        if(not proc in valuesDict):
          valuesDict[proc] = {}
  
        if(not type in valuesDict[proc]):
          valuesDict[proc][type] = []
  
        valuesDict[proc][type].append(float(time))
        haserror=0



#print "\n",

headers = []

it =  iter( valuesDict.keys())
for key in it:
  if(len(valuesDict[key].keys())>len(headers)):
    headers = valuesDict[key].keys()


print("%s, " % ("\"p\"".rjust(10,' ')), end='')
for header,last in lookahead(headers):
    if(not(last)):
      print("%s, " % (header.rjust(10,' ')), end='')
    else:
      print("%s" % (header.rjust(10,' ')))

for proc in sorted(valuesDict.keys()):
#  headers = 
  
  maxVals=[]
  for header,last in lookahead(valuesDict[proc].keys()):
    if(len(valuesDict[proc][header])>len(maxVals)):
      maxVals = valuesDict[proc][header]

  for i in range(0,len(maxVals)):  
    print("%s, " % (("%d" % int(proc)).rjust(10,' ')), end='')
  
    for type, last in lookahead(headers):
      val = float('nan')
      if(type in valuesDict[proc].keys()):
        if(i<len(valuesDict[proc][type])):
          val = valuesDict[proc][type][i]
      #print(val)
      if(not(last)):
        print("%s, " % (("%.5f" % val).rjust(10,' ')), end='')
      else:
        print("%s" % (("%.5f" % val).rjust(10,' ')))

    #print "%s %s %s %s" % ( str_proc.rjust(10,' ') , str_avg.rjust(10,' '), str_med.rjust(10,' '), str_stdev.rjust(10,' ') )

##compute the averages and stddev
#for proc in sorted(values.keys()):
#    val =  values[proc] 
#    avg = numpy.mean(val)
#    med = numpy.median(val)
#    stdev = numpy.std(val)
##    print "%d\t%10.5f\t%10.5f" % ( int(proc), avg, stdev )
#    str_proc = "%d" % int(proc)
#    str_avg = "%.5f" % avg
#    str_med = "%.5f" % med
#    str_stdev = "%.5f" % stdev
#    print("%s %s %s %s" % ( str_proc.rjust(10,' ') , str_avg.rjust(10,' '), str_med.rjust(10,' '), str_stdev.rjust(10,' ') ))
#



