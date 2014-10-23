#!/usr/bin/python

import sys
import re
import numpy

tregex = re.compile("^Factorization time: ([0-9.]+)")
pregex = re.compile("^[-]+([0-9]+)[-]+$")

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

values={}
with open(sys.argv[1], "rtU") as f:
  for line in f:
    m = pregex.match(line)
    if(m):
      proc = int(m.group(1))
    else:
      m = tregex.match(line)
      if(m):
        time = m.group(1)
        if(not proc in values):
          values[proc] = []

        values[proc].append(float(time))


#compute the averages and stddev
for proc in sorted(values.keys()):
    val =  values[proc] 
    avg = numpy.mean(val)
    stdev = numpy.std(val)
#    print "%d\t%10.5f\t%10.5f" % ( int(proc), avg, stdev )
    str_proc = "%d" % int(proc)
    str_avg = "%.5f" % avg
    str_stdev = "%.5f" % stdev
    print "%s %s %s" % ( str_proc.ljust(10,' ') , str_avg.ljust(10,' '), str_stdev.ljust(10,' ') )

