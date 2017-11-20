#!/bin/sh
# select particles from the known DF
smooth << eof
0      # for terminal input
#
1       # of mass components
n       # not a disk
plum    # halo type
1       # mass of halo
1       # radial scale of halo
6       # truncation radius
dejo    # DF type
0       # select an iso
200 50 20 # ne, nh, nhe
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv dejo200k.dfn run200.dfn
# create the .pcs0 file
mkpcs <<EOF
200
EOF
# set up the model ready to evolve
mpirun -np 5 pcs2dmp_mpi << eof
200
0
eof
# evolve the model
mpirun -np 5 galaxy_mpi <<EOF
200
0
EOF
