      real*8 function papotf( rs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Integrand for potential from a constant surface density annular ring of
c   material when the source point does not lie within the ring
c
c calling argument
      real*8 rs
c
c common block
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c local variables
      real*8 fr, fz, pot
c
      call rngfrc( rs, rf, zf, fr, fz, pot )
      papotf = pot * rs
      return
      end
