      real*8 function phgenf( rs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c regular integrand for phgend
c
c calling argument - the radius of the source
      real*8 rs
c
c common block
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c external
      real*8 gsigma
c
c local variables
      real*8 fr, fz, pot
c
      call rngfrc( rs, rf, zf, fr, fz, pot )
      phgenf = pot * rs * gsigma( rs )
      return
      end
