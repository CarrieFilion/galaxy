      subroutine qudini
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the coefficients for the Weinberg (1985) approximation to the
c   bar quadrupole field
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      integer idir
      common / ellipi / idir
c
      real*8 a1, a2, a3, axes( 3 ), ca( 3 )
      common / ellipr / a1, a2, a3, ca
      equivalence ( a1, axes( 1 ) )
c
c external
      external qudfun
      real*8 quad_inf
c
c local variables
      integer ier, inf
      real*8 bound, epsa, epsr
c
      if( .not. quapp )call crash( 'QUDINI',
     +                           'Quadrupole approximation not wanted' )
c set axes
      a1 = abar
      a2 = bbar
      a3 = cbar
c compute elliptic integrals
      do idir = 1, 2
        bound = 0
        inf = 1
        epsa = 1.e-10
        epsr = 1.e-8
        ca( idir ) = quad_inf( qudfun, bound, inf, epsa, epsr, ier )
        if( ier .gt. 0 )then
          print *, 'ier =', ier, ' from QUAD_INF'
          call crash( 'QUDINI', 'QUADPACK error' )
        end if
        ca( idir ) = a1 * a2 * a3 * ca( idir )
      end do
      alpha2 = 3 * a1**2 / ( 8. * a2 * a3 ) * ( ca( 2 ) - ca( 1 ) )
      beta2 = .4 * a1 * a2 * a3 * ( a2**2 - a1**2 ) /
     +                                             ( ca( 1 ) - ca( 2 ) )
      beta2 = beta2**0.2 / a1
      if( master )print *, 'alpha & beta', alpha2, beta2
      return
      end

      real*8 function qudfun( u )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
      real*8 u
c
c common blocks
c
      integer idir
      common / ellipi / idir
c
      real*8 a1, a2, a3, axes( 3 ), ca( 3 )
      common / ellipr / a1, a2, a3, ca
      equivalence ( a1, axes( 1 ) )
c
c local variables
      real*8 Delta
c
      Delta = ( a1**2 + u ) * ( a2**2 + u ) * ( a3**2 + u )
      Delta = sqrt( Delta )
      qudfun = 1. / ( Delta * ( axes( idir )**2 + u ) )
      return
      end
