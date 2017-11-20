      real*8 function gmtabd( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the mass interior to R for a disk model specified in
c    tabular form
c
c calling argument
      real*8 r
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
      integer i, iuse, lunit
      real*8 ar, rmin, rmax
      parameter ( lunit = 19 )
      save c, iuse, lam, rmin, rmax
c
      data iuse / 0 /
c
      if( iuse .ne. 1 )then
        call crash( 'GMTABD', 'Dummy version called' )
c open and read data file
        open( lunit, file = 'dsplin.dat', form = 'unformatted',
     +        status = 'old', iostat = i )
        if( i .ne. 0 )call crash( 'GMTABD', '.dat file not found' )
c read in data
        read( lunit )rad, mass
        close( lunit )
        rmin = rad( 2 )
        rmax = rad( m )
        ar = .5 * rmax
        gmtabd = splint2( rad, mass, m, ar, lam, c, .true. )
c reset flag
        iuse = 1
      end if
      ar = abs( r )
      if( ar .gt. rmax )then
        gmtabd = mass( m )
c interpolate
      else if( ar .ge. rmin )then
        gmtabd = splint2( rad, mass, m, ar, lam, c, .false. )
      else
c quadratic variation
        gmtabd = mass( 2 ) * ( ar / rad( 2 ) )**2
      end if
      return
      end
