c
c / model /
c
c maximum number of distinct mass components
      integer mcmp
      parameter ( mcmp = 10 )
c table sizes for isotropic DFs (King models and polytropes)
      integer nrdd
      parameter ( nrdd = 501 )
c maximum table size for adiabatic halo compression
      integer mradad
      parameter ( mradad = 10001 )
c
c / admin /
c
c maximum number of grids active at one time
      integer mgrids
      parameter ( mgrids = 4 )
c number of available solvers for gravitational field
      integer ncodes
      parameter ( ncodes = 11 )
c maximum number of nodes available for parallel processing
      integer mnodes
      parameter ( mnodes = 32 )
c maximum number of time-step zones
      integer mzones
      parameter ( mzones = 10 )
c maximum number of guard radii
      integer mguard
      parameter ( mguard = 10 )
c number of previous steps recorded for recentering
      integer npast
      parameter ( npast = 5 )
c pointer array for / mesh /
      integer mpt
      parameter ( mpt = 25 )
c maximum number of linked lists
      integer lists
      parameter ( lists = 257 )
c number of different analysis keywords
      integer ndattyp
      parameter ( ndattyp = 25 )
c
c / buffer /
c
c maximum number of particles in a group for time stepping
      integer mbuff
      parameter ( mbuff = 1000 )
c
c / grids /
c
c Fourier filter size
      integer lnfilt
      parameter ( lnfilt = 257 )
c maximum order of surface harmonic expansions
      integer s3maxl
      parameter ( s3maxl = 32 )
c
c / anlys /
c
c number of possible types of data in the .res file
      integer mtype
      parameter ( mtype = 28 )
c number of analys steps in one res file
      integer mtm
      parameter ( mtm = 10001 )
