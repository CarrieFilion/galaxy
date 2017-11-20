#!/bin/sh
# select particles from a known DF
smooth << eof
1      # input from script
100     # mass component
200 50 10 # numbers of particles
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv kaln100k.dfn run100.dfn
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
