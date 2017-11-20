      subroutine ajsetup
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Read values of roots and coefficients from appropriate data file
c   originally written by David Earn
c Modified by Jerry Sellwood and included in this package with permission
c
c common blocks containing basis choice and arrays of roots and coefficients
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
c local variables
      character ajfile*17, line*1000
      integer i, k, m, n
c
      if( basis .gt. mostkaj )call crash( 'AJSETUP',
     +                                             'mostkaj too small' )
c open the input file
      write( ajfile, '( a, i2.2, a )' )'../Basis/aj', basis, '.dat'
      open( 99, file = ajfile, status = 'old', form = 'formatted' )
c read the roots and coefficients for each required function
      do while ( .true. )
        read( 99, *, end = 1 )k, m, n
c check record header
        if( k .ne. basis )then
          if( master )print *, 'k =', k, '    basis =', basis
          call crash( 'AJSETUP', 'Wrong basis in file' )
        end if
c discard unneeded data
        if( ( m .gt. maxm ) .or. ( n .gt. maxn ) )then
          read( 99, '( a )' )line
          read( 99, '( a )' )line
        else
c check space
          if( m .gt. mostm )call crash( 'AJSETUP', 'mostm too small' )
          if( n .gt. mostn )call crash( 'AJSETUP', 'mostn too small' )
c read in data
          read( 99, '( a )' )line
          read( line, *, err = 2 )( root( m, n, i ), i = 1, n )
          read( 99, '( a )' )line
          read( line, *, err = 3 )( xcoeff( m, n, i ), i = 0, k )
        end if
      end do
c various error traps
    2 if( master )print 200, k, m, n
      call crash( 'AJSETUP', 'not enough root data' )
    3 if( master )print 200, k, m, n
      call crash( 'AJSETUP', 'not enough xcoeff data' )
c finish up
    1 close( 99 )
      if( master )print 200, k, m, n
c check that enough data were found
      if( ( m .lt. maxm ) .or. ( n .lt. maxn ) )call crash( 'AJSETUP',
     +                                     'not enough data in file' )
      return
  200 format( 'AJSETUP: last (k,m,n) was ( ',
     +        i2, ',', i2, ',', i3, ' )' )
      end
