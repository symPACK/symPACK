#!/usr/bin/python

import sys
import re
import numpy
import fileinput

tregex = re.compile("^Factorization time: ([0-9.]+)")
tregex2 = re.compile("^Factorization done in: ([0-9.]+)s")
pregex = re.compile("^[-]+([0-9]+)[-]+$")
pregex2 = re.compile("[ ]*NUMBER OF WORKING PROCESSES[ ]*=[ ]*([0-9]+)")

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

values={}
#with open(sys.argv[1], "rtU") as f:
#  for line in f:
for line in fileinput.input():
    m = pregex.match(line)
    m2 = pregex2.match(line)
    if(m):
      proc = int(m.group(1))
    elif(m2):
      proc = int(m2.group(1))
    else:
      m = tregex.match(line)
      m2 = tregex2.match(line)
      if(m):
        time = m.group(1)
        if(not proc in values):
          values[proc] = []

        values[proc].append(float(time))
      elif(m2):
        time = m2.group(1)
        if(not proc in values):
          values[proc] = []

        values[proc].append(float(time))



#compute the averages and stddev
for proc in sorted(values.keys()):
    val =  values[proc] 

    print key,val

#    avg = numpy.mean(val)
#    med = numpy.median(val)
#    stdev = numpy.std(val)
##    print "%d\t%10.5f\t%10.5f" % ( int(proc), avg, stdev )
#    str_proc = "%d" % int(proc)
#    str_avg = "%.5f" % avg
#    str_med = "%.5f" % med
#    str_stdev = "%.5f" % stdev
#    print "%s %s %s %s" % ( str_proc.ljust(10,' ') , str_avg.ljust(10,' '), str_med.ljust(10,' '), str_stdev.ljust(10,' ') )

