#!/bin/sh
# create snapshot files - run200.pcs0 already exists
for i in 25 50 75 100 ; do
   mpirun -np 5 dmp2pcs_mpi <<eof 
200  # run number
${i}  # time of the snaphot needed
    eof
done
# delete intermediate .dmp files as they have no futher use
for i in 0 25 50 75 ; do
    rm -f run200.*dmp${i}
done
# but the run200*dmp100 files could be needed to continue the simulation
rm -f run200.grt run200.flg
