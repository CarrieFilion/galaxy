#!/bin/sh
# determine the potential and density profiles of the bulge and halo
compress << eof
0
3
n
hern
0.3333
0.2
10
hern
0
n
isot
0.7
30
5.
eddi
y
kt
0.6667
1.
5.0
none
1.5
4.5
y
5
n
n
1 1
y
y
y
y
y
eof
ln -s compress.dat run200.cmp
# create the .dfn file for the bulge
smooth << eof
1      # 
200    # run no of .dat file
1      # select component 1
400 100 25 
y      # equal mass particles?
y      # OK to accept change of mass fraction
y
y
eof
mv comp1000k.dfn run200.df1
# create the .dfn file for the halo
smooth << eof
1      # 
200    # run no of .dat file
2      # select component 2
400 100 25 
y      # equal mass particles?
y      # OK to accept change of mass fraction
y
y
eof
mv comp1000k.dfn run200.df2
# setup the simulation
mkpcs << eof
200   # run number
eof
# create the .dmp files ready for galaxy
mpirun -np 5 pcs2dmp_mpi << eof
200   # run number
0     # time of .pcs file
eof
# run galaxy
mpirun -np 5 galaxy_mpi << eof 
200   # run number
0     # time of .dmp file
eof
