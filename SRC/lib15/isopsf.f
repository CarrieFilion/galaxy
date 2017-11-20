      subroutine isopsf( s, y, f )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c external routine need by ODE integrators for the run of psi in
c   isotropic spherical models such as polytropes and King models
c
c the 2nd order equation for Psi is split into
c
c    dPsi     X                   dX
c    ----- = ---     (1)   and    -- = -s**2 rho    (2)
c     ds     s**2                 ds
c
c with y( 1 ) = Psi and y( 2 ) = X = s**2 dPsi/ds
c
c calling arguments
      real*8 s, y( 2 ), f( 2 )
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external - density as function of psi
      real*8 isrhop
c
c local variable
      include 'inc/pi.f'
c
c gradient of equation (1)
      if( s**2 .gt. 0.d0 )then
        f( 1 ) = y( 2 ) / s**2
      else
        f( 1 ) = 0
      end if
      f( 2 ) = -s**2 * isrhop( y( 1 ) )
      return
      end
