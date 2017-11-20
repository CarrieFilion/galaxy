      real*8 function phigsz( r, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c gravitiational potential at any field point of a full mass Gaussian disc
c
c calling arguments
      real*8 r, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      common / gszphi / rad, zed
      real*8 rad, zed
c
c externals
      external phigzf
      real*8 quad_osc
c
c local variables
      integer ier
      real*8 b, epsa, epsr
c
      rad = r
      zed = abs( z )
c evaluate integral over many periods of Bessel function
      if( rad .gt. 1. )then
        b = 500. / sqrt( rad )
      else
        b = 500.
      end if
      epsa = 1.d-9
      epsr = epsa
c quadrature rule adapted for oscillatory integrands
      phigsz = quad_osc( phigzf, 0.d0, b, epsa, epsr, ier )
      if( ier .ne. 0 )then
        print *, 'ier =', ier, ' from DQAG'
        call crash( 'PHIGSZ', 'QUADPACK error' )
      end if
      return
      end
