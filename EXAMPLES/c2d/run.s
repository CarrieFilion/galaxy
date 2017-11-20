#!/bin/sh
smooth << eof
0
1            # 1 mass component
y            # this is a disk
mfk          # disk type
1            # mass scale
1            # radius scale
cpb0         # DF type
100 20 10    # Enter ne, nh and nhe
y            # equal mass particles
0            # random seed, or 0 for default seed
eof
mv cpb020k.dfn run100.dfn
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
