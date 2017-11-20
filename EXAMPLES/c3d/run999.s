#!/bin/sh
dfiter <<eof
999
1  # select cold start
n  # no, this is not a spherical model
y  # OK to overwrite old .pft file
eof
smooth << eof
1      # use file
999    # run #
100 20 10
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv dfit20k.dfn run999.dfn
# set up the model
mkpcs << eof 
999   # run number
y
eof
# create the .dmp files ready for galaxy
pcs2dmp << eof 
999   # run number
0     # time of .pcs file
eof
# run galaxy
galaxy << eof 
999   # run number
0     # time of .dmp file
eof
