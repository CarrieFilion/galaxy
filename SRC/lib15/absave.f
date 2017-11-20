      subroutine absave( m, nmodes, sumsq )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to store the best fit coefficients found by modfit into a
c   file called runXXX.mde for use by other programs (eg drwmode)
c
c calling arguments
      integer m, nmodes
      real sumsq
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
c local variables
      character*4 bstr( 5 )
      integer i, imode, istat, j, iuse
      real ahol
      save iuse
c
      data iuse / 0 /
      data bstr / 'LGSP', 'DANL', 'SFPC', 'SPHB', 'ZANL' /
c
c skip if file is already open
      if( iuse .eq. 0 )then
c try to open an old file
        j = -1
        call opnfil( j, 'mde', 'unformatted', 'old', 'append', istat )
        if( istat .ne. 0 )then
          print *, 'No old file found'
          call gtintg( 'Enter 0 to open a new file, 1 to stop', i )
          if( i .eq. 1 )stop
          print *, 'Opening new file'
          call opnfil( j, 'mde', 'unformatted', 'new', 'seq', istat )
c write header record
          call hedrec ( j, .false. )
        end if
        iuse = 1
      end if
c add extra record
      i = 0
      if( lgsp )i = 1
      if( danl )i = 2
      if( sfpc )i = 3
      if( sphb )i = 4
      if( zanl )i = 5
      if( i .eq. 0 )call crash( 'ABSAVE', 'Unrecognized data type' )
      read( bstr( i ), '( a4 )' )ahol
      write( j )irun, ahol, m, jt, kt, jp, kp, nmodes, apodc, sumsq,
     +           ncp
      write( j )( freq( 1, imode ), freq( 2, imode ),
     +             ( alp( i, imode ), bet( i, imode ), i = jp, kp ),
     +             imode = 1, nmodes )
      return
      end
