#!/bin/sh
# set up the model
mkpcs << eof 
100   # run number
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
