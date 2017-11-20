      subroutine sfpset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up the solution for the gravitational field using the SFP
c   approach.  It is adapted from David Earn's routine called SFPINPUT
c Modified by Jerry Sellwood and included in this package with permission
c
c This routine reads input for the Smooth-Field-Particle (SFP) code from the
c   input data file runxxx.dat, using subroutine getline.
c
c  common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local variables
      character eof*1, line*80, keyword*40
c
      if( .not. sf2d )call CRASH( 'sfpset', 'Wrong field method' )
c
c set the ascii end of file character for use with subroutine getline
      eof = char( 4 )
c name of the basis set to be used
      call getline( ni, line )
      if( line( 1:11 ) .eq. 'ABEL-JACOBI' )then
        basset = 'ablj'
      else if( line( 1:6 ) .eq. 'BESSEL' )then
        basset = 'bess'
      else if( line( 1:6 ) .eq. 'LOGSPI' )then
        basset = 'lgsp'
      else
        print*, 'Basis ', line( 1:40 ), ' is not available.'
        call crash( 'SFPSET', 'unknown basisname' )
      end if
c choose which of the Abel-Jacobi sets
      if( basset .eq. 'ablj' )then
        read( line( 12:14 ), '( i3 )' )basis
        if( basis .lt. 0 .or. basis .gt. mostkaj )then
          print*, 'You entered basis =', basis
          call crash( 'SFPSET', 'Basis out of range' )
        end if
      else if( basset .eq. 'bess' )then
        basis = 0
        read( line( 11:20 ), '( f10.4 )' )deltak
c softening length
        call getline( ni, line )
        read( line, '( 10x, f10.6 )' )softl
      else if( basset .eq. 'lgsp' )then
        basis = 0
        read( line( 11:20 ), '( f10.4 )' )deltak
      end if
c determine which values of n and m to use
      call getline( ni, line )
      if( line( 1:19 ) .ne. 'MINM MAXM MINN MAXN' )then
        print*, 'I do not understand ', line
        call crash( 'SFPSET', 'Unrecognized (m,n) input format' )
      end if
      call getline( ni, line )
      read( line, '( i4, 3(1x,i4) )' ) minm, maxm, minn, maxn
      if( ( minm .gt. maxm ) .or. ( minm .lt. 0 ) .or.
     +    ( minn .gt. maxn ) .or. ( minn .lt. 0 ) )then
        write( *, '( 3x, a )' )'minm maxm minn maxn'
        write( *, '( 2x, 4( 1x, i4 ) )' )minm, maxm, minn, maxn
        call crash( 'SFPSET', 'Min-max mixup' )
      else if( ( maxm .gt. mostm ) .or. ( maxn .gt. mostn ) )then
        print *, ' maxm =', maxm, '   maxn =', maxn
        print *, 'mostm =', mostm,'  mostn =', mostn
        print *, 'Increase these parameters and recompile the library'
        call crash( 'SFPSET', 'max m or n too large' )
      end if
c number of radii to tabulate functions
      call getline( ni, line )
      if( line( 1:5 ) .ne. 'NRTAB' )then
        print*, 'I do not understand ', line
        call crash( 'SFPSET', 'Not nrtab line' )
      end if
      read( line( 11:24 ), * )nrtab
      if( nrtab .le. 10 )then
        print *, 'value of nrtab read:', nrtab
        call crash( 'SFPSET', 'Nonsense value for nrtab' )
      end if
c Determine the minimum and maximum allowed value of maxr.  These are
c   entered as a fraction of the initial outer edge of the disc.
      if( basset .eq. 'ablj' )then
        call getline( ni, line )
        read( line, '( a7 )' )keyword
        if( keyword .ne. 'MINMAXR' )then
          print*, 'I do not understand ', line
          call crash( 'SFPSET', 'Bad minmaxr keyword' )
        end if
        read( line( 9:22 ), '( g14.6 )' )minmaxr
c
        call getline( ni, line )
        if( line .eq. eof )then
          print *, 'You must also choose a maximum value for MAXR.'
          print *, 'Enter it as a fraction of the initial outer ',
     +            'edge of the disc.'
          call crash( 'SFPSET', 'EOF before MAXMAXR specification' )
        end if
        read( line, '( a7 )' ) keyword
        if( keyword .ne. 'MAXMAXR' )then
          print *, 'I dont understand ', line
          call crash( 'SFPSET', 'bad maxmaxr keyword' )
        end if
        read( line( 9:22 ), '( g14.6 )' )maxmaxr
        maxmaxr = maxmaxr
c
        if( minmaxr .lt. 0. )call crash( 'SFPSET',
     +                                 'minmaxr is negative.' )
        if( minmaxr .eq. 0. )print *,
     +                  'maxr can have any value less than', maxmaxr
        if( maxmaxr .lt. minmaxr )call crash( 'SFPSET',
     +                                      'maxmaxr < minmaxr.' )
      else if( ( basset .eq. 'bess' ) .or. ( basset .eq. 'lgsp' ) )then
c set fiducial values
        minmaxr = 1.2
        maxmaxr = 1.2
      end if
      return
      end
