#!/bin/sh
rm -f run$1.lis
mpirun -np 4 ptest_mpi << eof 
$1
eof
cat run$1.lis
