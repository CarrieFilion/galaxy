#!/bin/sh
# select particles from a known DF
smooth << eof
0      # 
1      # 
n      # 
plum   # disk type
1      # mass
1      # radial scale
10     # rmax
dejo   # DF type 
0      # DF const
50 20 10
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv dejo10k.dfn run100.dfn
# create the .pcs file
mkpcs << eof 
100   # run number
eof
# create the .dmp files
/opt/openmpi/gnu/bin/mpirun -np 5 pcs2dmp_mpi << eof 
100   # run number
0     # start time
eof
# evolve the model
/opt/openmpi/gnu/bin/mpirun -np 5 galaxy_mpi << eof 
100   # run number
0     # start time
eof
