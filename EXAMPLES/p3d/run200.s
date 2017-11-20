#!/bin/sh
# select particles
smooth << eof
1     # read the .dat file
200   # run number
200 50 10
y     # equal mass particles
0     # use default seed
y     # OK to create new .dft file
y     # yes difference is acceptable
eof
mv shue100k.dfn run200.dfn
# set up a .pcs file
mkpcs << eof 
200   # run number
eof
# create the .dmp files ready for galaxy
pcs2dmp << eof 
200   # run number
0     # time of .pcs file
eof
# run galaxy
galaxy << eof 
200   # run number
0     # time of .dmp file
eof
