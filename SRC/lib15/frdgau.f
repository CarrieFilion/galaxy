      real*8 function frdgau( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c attraction in the mid-plane of a softened full-mass Gaussian disc
c   uses QUADPACK routine DQAG
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      common / gszphi / rad, zed
      real*8 rad, zed
c
c externals
      external frdgaf
      real*8 algrng2, quad_osc
c
c local arrays
      real*8, allocatable :: x( : )
      real*8, allocatable :: y( : )
c
c local variables
      integer i, ier, iuse, j, ntab
      real*8 b, epsa, epsr
      parameter ( ntab = 100 )
      save iuse, x, y
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
c create table
        print *, ' building a table for force from a Gaussian disc'
        zed = epsiln
        allocate ( x( ntab + 1 ) )
        allocate ( y( ntab + 1 ) )
        do j = 1, ntab
          i = j - 1
c uniform spacing in r / ( 1 + r )
          rad = dble( i ) / dble( ntab - i )
          x( j ) = dble( i ) / dble( ntab )
c evaluate integral over many periods of Bessel function
          if( rad .gt. 1. )then
            b = 500. / sqrt( rad )
          else
            b = 500.
          end if
          epsa = 1.d-9
          epsr = epsa
          y( j ) = quad_osc( frdgaf, 0.d0, b, epsa, epsr, ier )
          if( ier .ne. 0 )then
            print *, 'ier =', ier, ' from QUAD_OSC'
            call crash( 'FRDGAU', 'QUADPACK error' )
          end if
          y( j ) = y( j ) * ( 1. + rad )
        end do
c end value
        x( ntab + 1 ) = 1
        y( ntab + 1 ) = -1
        iuse = 1
        print *, 'table ready'
      end if
c look up value in table
      rad = abs( r )
      frdgau = algrng2( x, y, ntab + 1, rad / ( 1. + rad ) ) /
     +                                                      ( 1. + rad )
      frdgau = sign( frdgau, -r )
      return
      end
