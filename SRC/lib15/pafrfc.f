      real*8 function pafrfc( rs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Integrand for Cauchy principal value integral for radial force from a
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
      real*8 pafrfn
c
      pafrfc = ( rs - rf ) * pafrfn( rs )
      return
      end
