#!/bin/sh
# find the DF that is in equilibrium when a massive disk is added
compress <<EOF
0
2
# component 1
y              # a disc
exp            # type keyword
1              # mass in this pop
1              # length scale of this pop
5              # truncation radius of this pop
none           # no DF
1.5            # initial Q
4.5            # inner edge of outer taper
# component 2
n              # not a disk
king           # halo type
3              # halo parameter
8              # halo mass       
10             # halo length scale
king           # halo DF type
#
y              # accept disk vertical profile
0              # suppress graphics output
y              # OK to overwite (needed only if there was a false start)
EOF
ln -s compress.dat run100.cmp
# select particles from this DF
smooth << eof
1      # read data form .dat file
100    # run number
500 100 50 # ne, nh, nhe
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv comp2500k.dfn run100.dfn
# create the .pcs file
mkpcs << eof 
100   # run number
eof
# create the .dmp files ready for galaxy
mpirun -np 5 pcs2dmp_mpi << eof
100   # run number
0     # time of .pcs file
eof
# run galaxy
mpirun -np 5 galaxy_mpi << eof 
100   # run number
0     # time of .dmp file
eof
