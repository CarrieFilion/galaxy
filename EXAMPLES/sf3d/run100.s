#!/bin/sh
# select particles from a known DF
smooth << eof
1      # input from .dat file
100    # run number of .dat file
100 20 5
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv kaln10k.dfn run100.dfn
# create the .pcs file
mkpcs << eof 
100   # run number
eof
# create the .dmp files
pcs2dmp << eof 
100   # run number
0     # start time
eof
# evolve the model
galaxy << eof 
100   # run number
0     # start time
eof
