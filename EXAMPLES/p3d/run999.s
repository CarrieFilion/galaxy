#!/bin/sh
# find an equilibrium DF for the bulge
dfiter <<eof
999
1  # select cold start
n  # this is not a spherical model
y  # OK to overwrite old .pft file
eof
# select particles from the equilibriumn DF
smooth << eof
1
999         # run number
100 50 20   # need 100K particles
y           # equal mass
0           # accept default random seed
eof
mv dfit100k.dfn run999.dfn
# set up a .pcs file
mkpcs << eof 
999   # run number
y     # yes, accept change
eof
# create the .dmp files ready for galaxy
mpirun -np 5 pcs2dmp_mpi << eof 
999   # run number
0     # time of .pcs file
eof
# run galaxy
mpirun -np 5 galaxy_mpi << eof 
999   # run number
0     # time of .dmp file
eof
