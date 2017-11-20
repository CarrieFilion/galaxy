#!/bin/sh
# select particles from the DF
smooth <<eof
0     # input from script
1     # 1 mass component
n     # not a disk
unis  # mass type
1     # mass of this component
1     # radial scale of this component
usps  # DF type
200 50 20 # generate 200K particles
y     # equal mass particles
0     # accept default random seed
eof
mv usps200k.dfn run100.dfn
# set up a .pcs file
mkpcs << eof 
100   # run number
eof
# create the .dmp files ready for galaxy
pcs2dmp << eof
100
0
eof
# evolve the model
galaxy <<eof
100
0
eof
