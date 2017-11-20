      subroutine rholm( r, wt )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the surface harmonic coefficients for a given density distribution
c   on a shell of radius r (in natural units)
c
c calling arguments
      real r, wt( 2, * )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      common / lmrho / l, m, r2, theta2
      integer l, m
      real*8 r2, theta2
c
c externals
      real*8 rhof1i, rhof1r
c
c local arrays
      real*8, allocatable :: absc( : ), wgt( : )
c
c local variables
      integer i, k, npts
      real*8 x, zero
      parameter ( npts = 24, zero = 0 )
      include 'inc/pi.f'
c
      r2 = r
c      rho0 = dskrho( r2, zero, zero )
      k = 0
      do l = 0, s3lmax
        do m = 0, l
          k = k + 1
          wt( 1, k ) = 0
          wt( 2, k ) = 0
c          if( rho0 .gt. zero )then
c odd l terms are zero for mass distributions symmetric about the plane
            if( mod( l, 2 ) .eq. 0 )then
              if( .not. allocated( absc ) )then
                allocate ( absc( npts ) )
                allocate ( wgt( npts ) )
              end if
              call GLqtab( zero, .5 * pi, npts, wgt, absc )
              wt( 1, k ) = 0
              do i = 1, npts
                x = absc( i )
                wt( 1, k ) = wt( 1, k ) + wgt( i ) * rhof1r( x )
              end do
              wt( 1, k ) = 2. * wt( 1, k )
              if( m .gt. 0 )then
                wt( 2, k ) = 0
                do i = 1, npts
                  x = absc( i )
                  wt( 2, k ) = wt( 2, k ) + wgt( i ) * rhof1i( x )
                end do
                wt( 2, k ) = 2. * wt( 2, k )
              end if
            end if
c          end if
        end do
      end do
      return
      end

      real*8 function rhof1r( theta )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for outer integral (over polar angle) for computation of real rholm
c   uses QUADPACK routine DQAG
c
c calling argument
      real*8 theta
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      common / lmrho / l, m, r2, theta2
      integer l, m
      real*8 r2, theta2
c
c externals
      external rhof2r
      real*8 quad_osc
c
c local variables
      integer ier
      real*8 epsa, epsr, zero
      parameter ( zero = 0 )
      include 'inc/pi.f'
c
c preserve calling argument
      theta2 = theta
      call s3tplm( cos( theta2 ), .false. )
c
      epsa = 1.d-5
      epsr = epsa
      rhof1r = quad_osc( rhof2r, zero, 2. * pi, epsa, epsr, ier )
      if( ier .ne. 0 )then
        print *, 'ier =', ier, ' from DQAG'
        call crash( 'RHOF1R', 'QUADPACK error' )
      end if
      rhof1r = rhof1r * sin( theta2 )
      return
      end

      real*8 function rhof2r( phi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for inner integral (over azimuth) for computation of real rholm
c
c calling argument
      real*8 phi
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      common / lmrho / l, m, r2, theta2
      integer l, m
      real*8 r2, theta2
c
c external
      real*8 dskrho
c
c local variables
      integer k
      real*8 arg, rad, x, y, z
c
      z = r2 * cos( theta2 )
      rad = r2 * sin( theta2 )
      x = rad * cos( phi )
      y = rad * sin( phi )
      k = ( l * ( l + 1 ) ) / 2 + m + 1
      arg = m * phi
      rhof2r = plm( k ) * cos( arg ) * dskrho( x, y, z )
      return
      end

      real*8 function rhof1i( theta )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for outer integral (over polar angle) for computation of
c   imaginary rholm
c   uses QUADPACK routine DQAG
c
c calling argument
      real*8 theta
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      common / lmrho / l, m, r2, theta2
      integer l, m
      real*8 r2, theta2
c
c externals
      external rhof2i
      real*8 quad_osc
c
c local variables
      integer ier
      real*8 epsa, epsr, zero
      parameter ( zero = 0 )
      include 'inc/pi.f'
c
c preserve calling argument
      theta2 = theta
      call s3tplm( cos( theta2 ), .false. )
c
      epsa = 1.d-5
      epsr = epsa
      rhof1i = quad_osc( rhof2i, zero, 2. * pi, epsa, epsr, ier )
      if( ier .ne. 0 )then
        print *, 'ier =', ier, ' from DQAG'
        call crash( 'RHOF1I', 'QUADPACK error' )
      end if
      rhof1i = rhof1i * sin( theta2 )
      return
      end

      real*8 function rhof2i( phi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for inner integral (over azimuth) for computation of
c   imaginary rholm
c
c calling argument
      real*8 phi
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      common / lmrho / l, m, r2, theta2
      integer l, m
      real*8 r2, theta2
c
c external
      real*8 dskrho
c
c local variables
      integer k
      real*8 arg, rad, x, y, z
c
      z = r2 * cos( theta2 )
      rad = r2 * sin( theta2 )
      x = rad * cos( phi )
      y = rad * sin( phi )
      k = ( l * ( l + 1 ) ) / 2 + m + 1
      arg = m * phi
      rhof2i = plm( k ) * sin( arg ) * dskrho( x, y, z )
      return
      end
