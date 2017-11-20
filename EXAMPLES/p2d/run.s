#!/bin/sh
# set up a smooth particle disk with circular orbits
cold << eof
0
1
y
soft
0.367879
1
1
8
kaln
999
0
6
5000
y
0
eof
mv kaln5k.dfn run100.dfn
# create the .pcs file
mkpcs << eof 
100   # run number
eof
# create the .dmp file(s) ready for galaxy
pcs2dmp << eof
100   # run number
0     # time of .pcs file
eof
# run galaxy
galaxy << eof
100   # run number
0     # time of .dmp file
eof
