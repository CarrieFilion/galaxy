      real*8 function frdisi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to evaulate the central attraction at radius r in the plane of
c   a thin disc with a surface density distribution given by gsigmt
c   (gsigmt is an accelerated version of gsigmi, which is the disc surface
c   obtained from a double integral over the DF)
c   uses QUADPACK routines DQAGS and DQAWC
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
      common / feildp / rf, zf
      real*8 rf, zf
c
c externals - frdsfn is the basic integrand, frdsfc is needed for the
c   Cauchy principal value part
      external frdsfc, frdsfn
      real*8 algrng2, quad_chy, quad_gnr
c
c local arrays
      integer ntab
      parameter ( ntab = 100 )
      real*8, allocatable :: x( : )
      real*8, allocatable :: y( : )
c
c local variables
      integer i, ier, iuse, j
      real*8 epsa, epsr, rad
      save iuse, x, y
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
c create table
        print *, ' building a table of radial force from the disc'
        zf = epsiln
        allocate ( x( ntab + 1 ) )
        allocate ( y( ntab + 1 ) )
        do j = 1, ntab
          i = j - 1
c uniform spacing in r / ( 1 + r )
          rad = dble( i ) / dble( ntab - i )
          x( j ) = dble( i ) / dble( ntab )
          if( rad .gt. 0.d0 )then
c integrand is regular everywhere for a field point outside the disc
            if( ( rad .lt. rhole - 1.d-8 ) .or.
     +          ( rad .gt. rmax + 1.d-8 ) .or.
     +          ( epsiln .gt. 0.d0 ) )then
              rf = rad
              epsa = 1.d-5
              epsr = epsa
              y( j ) = quad_gnr( frdsfn, rhole, rmax, epsa, epsr, ier)
              if( ier .gt. 0 )then
                print *, 'ier =', ier, ' from QUAD_GNR'
                call crash( 'FRDISI', 'QUADPACK error' )
              end if
            else
c Cauchy principal value needed when field point is within the disc
              epsa = 1.d-5
              epsr = epsa
              ier = 0
              rf = min( rad, rmax - 1.d-8 )
              rf = max( rf, rhole + 1.d-8 )
              y( j ) =
     +              quad_chy( frdsfc, rhole, rmax, rf, epsa, epsr, ier )
              if( ier .gt. 0 )then
                print *, 'ier =', ier, ' from QUAD_CHY'
                call crash( 'FRDISI', 'QUADPACK error' )
              end if
            end if
            y( j ) = y( j ) * ( 1. + rad )
          end if
        end do
c end values
        y( 1 ) = 0
        x( ntab + 1 ) = 1
        y( ntab + 1 ) = -1
        iuse = 1
        print *, 'table ready'
      end if
c look up value in table
      rad = abs( r )
      frdisi = algrng2( x, y, ntab + 1, rad / ( 1. + rad ) ) /
     +                                                      ( 1. + rad )
      if( r .lt. 0.d0 )frdisi = -frdisi
      return
      end
