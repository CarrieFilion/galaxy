      subroutine sf3set
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up the solution for the gravitational field using the SFP
c   approach.  It reads input from the data file runXXX.dat.
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
      character eof*1, line*80
c
      if( .not. sf3d )call CRASH( 'SF3SET', 'Wrong method' )
c
c set the ascii end of file character for use with subroutine getline
      eof = char( 4 )
c name of the basis set to be used
      call getline( ni, line )
      if( line( 1:6 ) .eq. 'BESSEL' )then
        basset = 'bess'
      else
        print *, 'Basis ', line( 1:40 ), ' is not available.'
        call crash( 'SF3SET', 'Unknown basis type' )
      end if
c read parameters
      if( basset .eq. 'bess' )then
        read( line( 11:20 ), '( f10.4 )' )deltak
c softening length
        call getline( ni, line )
        read( line, '( 10x, f10.6 )' )softl
      end if
c determine which values of n and m to use
      call getline( ni, line )
      if( line( 1:19 ) .ne. 'MINM MAXM MINN MAXN' )then
        print*, 'I do not understand ', line
        call crash( 'SF3SET', 'Unknown (m,n) input format' )
      end if
      call getline( ni, line )
      read( line, '( i4, 3(1x,i4) )' ) minm, maxm, minn, maxn
      if( ( minm .gt. maxm ) .or. ( minm .lt. 0 ) .or.
     +    ( minn .gt. maxn ) .or. ( minn .lt. 0 ) )then
        write( *, '( 3x, a )' )'minm maxm minn maxn'
        write( *, '( 2x, 4( 1x, i4 ) )' )minm, maxm, minn, maxn
        call crash( 'SF3SET', 'Min-max mixup' )
      else if( ( maxm .gt. mostm ) .or. ( maxn .gt. mostn ) )then
        print*, ' maxm =', maxm, '   maxn =', maxn
        print*, 'mostm =', mostm,'  mostn =', mostn
        call crash( 'SF3SET', 'Max m or n too large' )
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
        call crash( 'SF3SET', 'Nonsense value for nrtab' )
      end if
c set fiducial values
      minmaxr = 1.2
      maxmaxr = 1.2
c number of planes
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'PLAN' )then
        print *, 'next card is'
        print *, line
        call crash( 'SF3SET', 'nplanes card not found' )
      end if
      read( line, '( 10x, i10 )' )ngz
      read( line, '( 20x, f10.4 )' )dzg
      return
      end
