#!/bin/sh
# run compress to create the equilibrium halo DF
compress << eof
0       # for terminal input
2       # number of components
#
n       # it is not a disk
plum    # component type
1       # mass
1       # scale length
10      # truncation radius
dejo    # DF type
0       # DF const (isotropic)
#
n       # it is not a disk
hern    # component type
0.2     # mass of this component
0.1     # radial scale
15      # truncation radius
none    # DF type
#
0       # no graphics output
y       # OK to overwrite
eof
ln -s compress.dat run300.cmp
# select particles from the halo DF
smooth << eof
1       # for file input
300     # run number of .dat file
200 50 20 # ne, nh, nhe
y       # equal mass particles?
0       # use default random seed
y       # OK to accept change of mass fraction
eof
mv comp200k.dfn run300.dfn
# create the .pcs0 file
mkpcs <<eof
300
eof
# set up the model ready to evolve
mpirun -np 5 pcs2dmp_mpi << eof
300
0
eof
# evolve the model
mpirun -np 5 galaxy_mpi <<eof
300
0
eof
