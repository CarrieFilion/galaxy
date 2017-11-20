      real*8 function actj1( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the radial action for a planar orbit
c
c calling arguments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      include 'inc/orbval.f'
c
c externals
      real*8 akappa, Lzmax, Phitot, rcirc
c
c local arrays
      integer npts
      parameter ( npts = 64 )
      real*8 wgt( npts ), acs( npts )
c
c local variables
      integer i
      logical set
      real*8 a, r
      include 'inc/pi.f'
c
c check arguments
      if( Lz - Lzmax( E ) .gt. 1.d-5 )then
        print *, 'E, Lz  = ', E, Lz
        call crash( 'ACTJ1', 'Impossible arguments' )
      end if
      set = .false.
      actj1 = 0
      if( ncmp .eq. 1 )then
c use analytic expression where available
        if( ctype( 1 ) .eq. 'ISOC' )then
c isochrone disc
          a = 0.
          actj1 = 1. / sqrt( -2. * E ) - .5 * abs( Lz ) -
     +                                        sqrt( .25 * Lz * Lz + 1. )
          actj1 = max( actj1, a )
          set = .true.
        end if
      end if
      if( .not. set )then
c find peri and apo
        call rlims( E, Lz )
        a = rapo - rperi
c check for almost circular orbits
        if( a .gt. 0. )then
          if( a / ( rperi + rapo ) .lt. 1.d-5 )then
            r = rcirc( Lz )
            a = .5 * a
            actj1 = .5 * a**2 * akappa( r )
          else
c defining expression - Gauss-Legendre quadradure from peri to apo
            call GLqtab( rperi, rapo, npts, wgt, acs )
c form integral round half radial oscillation
            actj1 = 0
            do i = 1, npts
              r = acs( i )
              a = 2. * ( E - Phitot( r ) ) - ( Lz / r )**2
              if( a .gt. 0.d0 )actj1 = actj1 + wgt( i ) * sqrt( a )
            end do
c normalize - integral was over half the radial orbit
            actj1 = actj1 / pi
          end if
        end if
      end if
      return
      end
