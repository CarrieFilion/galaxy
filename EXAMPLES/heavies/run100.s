#!/bin/sh
# select 20K particles from the DF of an isotropic Plummer sphere
smooth << eof
0      # input from script
1      # number of mass components
n      # it is not a disk
plum   # halo type
1      # mass
1      # radial scale
6      # truncation radius
dejo   # DF type
0      # no anisotropy
100 20 10
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv dejo20k.dfn run100.df1
# create a 2nd file with 100 particles selected from the same DF
smooth << eof
0      # input from script
1      # number of mass components
n      # it is not a disk
plum   # halo type
1      # mass
1      # radial scale
6      # truncation radius
dejo   # DF type
0      # no anisotropy
10 5 2
y      # equal mass particles?
0      # use default random seed
y      # OK to accept change of mass fraction
eof
mv dejo0k.dfn run100.df2
# create the .pcs file
mkpcs << eof 
100   # run number
n     # do not rescale velocities
y     # yes, this is acceptable
n     # do not rescale velocities
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
