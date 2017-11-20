      subroutine scfset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up the solution for the gravitational field using the SCF
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
      if( .not. scf )call crash( 'SCFSET', 'Wrong method' )
c
c set the ascii end of file character for use with subroutine getline
      eof = char( 4 )
c name of the basis set to be used
      call getline( ni, line )
      if( line( 1:6 ) .eq. 'HERNQU' )then
        basset = 'hern'
      else if( line( 1:6 ) .eq. 'PLUMME' )then
        basset = 'plum'
      else
        print *, 'Basis ', line( 1:40 ), ' is not available.'
        call crash( 'SCFSET', 'Unknown basis type' )
      end if
c radius of volume for analysis and plotting etc
      read( line( 11:80 ), *, err = 1 )rgrid( jgrid )
c choose azimuthal order
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'LMAX' )then
        print *, 'Next card is'
        print *, line
        call crash( 'SCFSET', 'lmax card not found' )
      end if
      read( line( 11:80 ), *, err = 1 )s3lmax
c choose radial order
      call getline( ni, line )
      if( line( 1:4 ) .ne. 'NMAX' )then
        print *, 'Next card is'
        print *, line
        call crash( 'SCFSET', 'nmax card not found' )
      end if
      read( line( 11:80 ), *, err = 1 )nbas
      return
    1 print '( a )', line( 1:60 )
      call crash( 'SCFSET', 'Error reading line' )
      stop
      end
