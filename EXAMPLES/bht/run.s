#!/bin/sh
# set up 10k particles from an isotropic Plummer sphere
smooth << eof
0        # for terminal input
#
1        # of mass components
n        # not a disk
plum     # halo type
1        # mass of halo
1        # radial scale of halo
5        # truncation radius
dejo     # DF type
0        # isotropic DF
100 20 5 # ne, nh, nhe
y        # equal mass particles?
0        # use default random seed
y        # OK to accept change of mass fraction
eof
mv dejo10k.dfn run100.dfn
# create the .pcs file
mkpcs << eof 
100   # run number
eof
# create the .dmp files ready for galaxy
pcs2dmp << eof 
100   # run number
0     # time of .pcs file
eof
# run galaxy
galaxy << eof 
100   # run number
0     # time of .dmp file
eof
