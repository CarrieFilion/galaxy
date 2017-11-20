c radial arrays for supplementary forces (to compensate for softening etc.)
      integer mrsup, nrsup
      parameter ( mrsup = 502 )
      logical lsupst
      real alp2, hrfac, htab( 3, mrsup ), selfe
      common / frcsup / nrsup, lsupst, hrfac, selfe, htab, alp2
c
      real hpot( mrsup )
      common / pothal / hpot
c
      real self( 12, mrsup )
      common / selfe / self
c
      real sff( 4, mrsup )
      common / slffor / sff
