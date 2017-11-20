      real*8 function papotc( rs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Integrand for Cauchy principal value integral for potential of a
c   constant surface density annular ring of mass.
c
c calling argument
      real*8 rs
c
c common block
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c external - integrand without the ( rs - rf ) term
      real*8 papotf
c
      papotc = ( rs - rf ) * papotf( rs )
      return
      end
