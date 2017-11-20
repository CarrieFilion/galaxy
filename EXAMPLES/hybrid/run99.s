#!/bin/sh
# create the .pcs file
mkpcs << eof 
99    # run number
eof
# create the .dmp files ready for galaxy
mpirun -np 5 pcs2dmp_mpi << eof 
99    # run number
0     # time of .pcs file
eof
# run galaxy
mpirun -np 5 galaxy_mpi << eof 
99    # run number
0     # time of .dmp file
eof
