#!/bin/sh
# set up a smooth particle disk with circular orbits
smooth << eof
0
1   # 1 mass component
y    # 'tis a disk
mtz   # disk type
1     # mass scale
1     # length scale
6     # trunk radius
zang    #DF type
1.5 # initial Q
6   # inner edge outer taper?
100 100 60 # ne, nh, nhe (number of particles?)
y  # equal mass particles
0      #  random seed, 0 for default seed
y      # ok to accept change of mass fn ?
eof
mv zang600k.dfn run400.dfn
# create the .pcs file
mkpcs << eof
400   # run number
eof
# create the .dmp file(s) ready for galaxy
pcs2dmp << eof
 400   # run number
 0     # time of .pcs file
 eof
# run galaxy
galaxy << eof
 400   # run number
 0     # time of .dmp file
eof

