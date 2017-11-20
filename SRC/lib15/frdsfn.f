      real*8 function frdsfn( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to evaluate the radial force contribution at the field point
c  stored in / fieldp / from an axisymmetric mass element at r
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      common / feildp / rf, zf
      real*8 rf, zf
c
      include 'inc/model.f'
c
c externals - gsigmt is an accelerated version of gsigmi
      real*8 gsigma, gsigmt
c
      real*8 fr, fz, pot
      include 'inc/pi.f'
c
c routine to determine force components and potential
      call rngfrc( r, rf, zf, fr, fz, pot )
c other parts of integrand
      if( ctype( icmp ) .eq. 'SPLY' )then
        frdsfn = 2. * pi * r * gsigma( r ) * fr
      else
        frdsfn = 2. * pi * r * gsigmt( r ) * fr
      end if
      return
      end
