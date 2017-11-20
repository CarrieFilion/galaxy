      real*8 function gmtabh( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the mass interior to r for a spherical model specified in
c    tabular form
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 splint2
c
c local arrays
      integer m
      parameter ( m = 1001 )
      real*8 rad( m ), mass( m ), lam( m + 4 ), c( m + 4 )
c
c local variables
      integer i, iuse, lunit, n
      real*8 ar, rmn, rmx
      parameter ( lunit = 19 )
      save c, iuse, lam, n, rmn, rmx
c
      data iuse / 0 /
c
      if( iuse .ne. 1 )then
        call crash( 'GMTABH', 'Dummy version called' )
c open and read data file
        open( lunit, file = 'hsplin.dat', form = 'unformatted',
     +        status = 'old', iostat = i )
        if( i .ne. 0 )call crash( 'GMTABH', '.dat file not found' )
c read in data
        read( lunit )rad, mass
        close( lunit )
        rmn = rad( 1 )
        rmx = rad( m )
        n = 0
        do i = 1, m, 10
          n = n + 1
          rad( n ) = rad( i )
          mass( n ) = cmpmas( 1 ) * mass( i )
        end do
        ar = .5 * rmx
        gmtabh = splint2( rad, mass, n, ar, lam, c, .true. )
c reset flag
        iuse = 1
      end if
      ar = abs( r )
      if( ar .gt. rmx )then
c outside tabulated range
        gmtabh = mass( n )
      else if( ar .ge. rad( 2 ) )then
c interpolate
        gmtabh = splint2( rad, mass, n, ar, lam, c, .false. )
      else
c quadratic variation - assumes a r^{-1} cusp
        gmtabh = mass( 2 ) * ( ar / rad( 2 ) )**2
      end if
      return
      end
