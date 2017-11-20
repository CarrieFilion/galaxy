# Copyright (C) 2015, Jerry Sellwood
SHELL = /bin/sh
#
# Fortran compiler for single processor code - may need to be changed
f90 = gfortran -w
#
# Fortran compiler for MPI code - may need to be changed
mf90 = mpif90
#
# edit this line to give the actual location of the PGPLOT library
#    PGPLOT may also require an X11 library
pgplot_lib = /usr/local/lib/pgplot/libpgplot.a /usr/X11/lib/libX11.dylib

#
# a bunch of system libraries that are needed by the linker
#    probably needs revision
sys_lib = -lgfortran -lpng12 -L/usr/X11R6/lib64 -L/usr/lib64 -lX11 
#
# The remainder of the file SHOULD NOT BE CHANGED
#
basic_lib = ../../lib/main15.a ../../lib/RAJ.a ../../lib/libutils.a
mypi_lib = ../../lib/main15_mpi.a
#
alls_lib = $(basic_lib) $(pgplot_lib) $(sys_lib)
allp_lib = $(mypi_lib) $(basic_lib) $(sys_lib)
analys cold compress corrplt dfiter dflook dmp2pcs estfreq gadget2pcs galaxy \
       genplt isotropy merge mkpcs modefit pcs2dmp plotpft ptest smooth weed \
       tipsy2pcs dmp2pcs_mpi galaxy_mpi pcs2dmp_mpi ptest_mpi : \
       aarrays.o analys.o cold.o compress.o corrplt.o dfiter.o dflook.o \
       dmp2pcs.o estfreq.o galaxy.o genplt.o gadget2pcs.o isotropy.o merge.o \
       mkpcs.o modefit.o pcs2dmp.o plotpft.o ptest.o smooth.o weed.o \
       inc tipsy2pcs.o pcs2dmp_mpi.o dmp2pcs_mpi.o galaxy_mpi.o ptest_mpi.o
	$(f90) analys.o $(alls_lib) -o analys
	$(f90) cold.o $(basic_lib) -o cold
	$(f90) compress.o $(alls_lib) -o compress
	$(f90) corrplt.o $(alls_lib) -o corrplt
	$(f90) dfiter.o $(basic_lib) -o dfiter
	$(f90) dflook.o $(alls_lib) -o dflook
	$(f90) dmp2pcs.o $(basic_lib) -o dmp2pcs
	$(f90) estfreq.o $(alls_lib) -o estfreq
	$(f90) gadget2pcs.o $(basic_lib) -o gadget2pcs
	$(f90) galaxy.o $(basic_lib) -o galaxy
	$(f90) genplt.o $(alls_lib) -o genplt
	$(f90) isotropy.o $(alls_lib) -o isotropy
	$(f90) merge.o $(basic_lib) -o merge
	$(f90) mkpcs.o $(basic_lib) -o mkpcs
	$(f90) modefit.o $(alls_lib) -o modefit
	$(f90) pcs2dmp.o $(basic_lib) -o pcs2dmp
	$(f90) plotpft.o $(alls_lib) -o plotpft
	$(f90) ptest.o $(basic_lib) -o ptest
	$(f90) smooth.o $(basic_lib) -o smooth
	$(f90) tipsy2pcs.o $(alls_lib) -o tipsy2pcs
	$(f90) weed.o $(basic_lib) -o weed
	$(mf90) dmp2pcs_mpi.o $(allp_lib) -o dmp2pcs_mpi
	$(mf90) galaxy_mpi.o $(allp_lib) -o galaxy_mpi
	$(mf90) pcs2dmp_mpi.o $(allp_lib) -o pcs2dmp_mpi
	$(mf90) ptest_mpi.o $(allp_lib) -o ptest_mpi
inc :
	ln -s ../../SRC/lib15/inc
aarrays.o : aarrays.f
	$(f90) -c aarrays.f
analys.o : analys.f inc
	$(f90) -c analys.f
cold.o : cold.f inc
	$(f90) -c cold.f
compress.o : compress.f inc
	$(f90) -c compress.f
corrplt.o : corrplt.f inc
	$(f90) -c corrplt.f
dfiter.o : dfiter.f inc
	$(f90) -c dfiter.f
dflook.o : dflook.f inc
	$(f90) -c dflook.f
dmp2pcs.o : dmp2pcs.f inc
	$(f90) -c dmp2pcs.f
estfreq.o : estfreq.f inc
	$(f90) -c estfreq.f
gadget2pcs.o : gadget2pcs.f inc
	$(f90) -c gadget2pcs.f
galaxy.o : galaxy.f inc
	$(f90) -c galaxy.f
genplt.o : genplt.f inc
	$(f90) -c genplt.f
isotropy.o : isotropy.f inc
	$(f90) -c isotropy.f
merge.o : merge.f inc
	$(f90) -c merge.f
mkpcs.o : mkpcs.f inc
	$(f90) -c mkpcs.f
modefit.o : modefit.f inc
	$(f90) -c modefit.f
pcs2dmp.o : pcs2dmp.f inc
	$(f90) -c pcs2dmp.f
plotpft.o : plotpft.f inc
	$(f90) -c plotpft.f
ptest.o : ptest.f inc
	$(f90) -c ptest.f
smooth.o : smooth.f inc
	$(f90) -c smooth.f
tipsy2pcs.o : tipsy2pcs.f inc
	$(f90) -c tipsy2pcs.f
weed.o : weed.f inc
	$(f90) -c weed.f
dmp2pcs_mpi.o : dmp2pcs_mpi.f inc
	$(mf90) -c aarrays.f dmp2pcs_mpi.f
galaxy_mpi.o : galaxy_mpi.f inc
	$(mf90) -c aarrays.f galaxy_mpi.f
pcs2dmp_mpi.o : pcs2dmp_mpi.f inc
	$(mf90) -c aarrays.f pcs2dmp_mpi.f
ptest_mpi.o : ptest_mpi.f inc
	$(mf90) -c aarrays.f ptest_mpi.f
clean :
	mv analys cold compress corrplt dfiter dflook dmp2pcs estfreq \
         gadget2pcs galaxy genplt isotropy merge mkpcs modefit pcs2dmp plotpft \
         ptest smooth tipsy2pcs weed \
         dmp2pcs_mpi galaxy_mpi pcs2dmp_mpi ptest_mpi ../../bin
	rm analys.o cold.o compress.o corrplt.o dfiter.o dflook.o dmp2pcs.o \
         estfreq.o gadget2pcs.o galaxy.o genplt.o isotropy.o merge.o mkpcs.o \
         modefit.o pcs2dmp.o plotpft.o ptest.o smooth.o weed.o \
         tipsy2pcs.o aarrays.o aarrays.mod mesh_data.mod inc  \
         dmp2pcs_mpi.o galaxy_mpi.o pcs2dmp_mpi.o ptest_mpi.o
#
# The author is very inexperienced with Makefiles and would welcome any
#   advice that would facilitate easier use of this make 

