      subroutine getset
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c general start up routine that sets run number and field method type:
c   sets flags for single or parallel processing
c   reads in the run number
c   verifies that the appropriate .dat file is present and opens it
c   verifies that .dat file contains the input run number
c   reads method type from the .dat file and calls codeset
c   opens .lis file with access = 'append' if possible
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
      integer iunit( 6 + 2 * mgrids )
      equivalence ( ni, iunit( 1 ) )
c
      include 'mpif.h'
c
c local variables
      character line*80
      character*3 code
      integer i, ierr, j
c
c initialize parallel processing
      call mpi_init( ierr )
      call mpi_comm_rank( mpi_comm_world, myid, ierr )
      call mpi_comm_size( mpi_comm_world, numprocs, ierr )
c set flags for parallel operation
      master = ( myid .eq. 0 ) .or. ( numprocs .eq. 1 )
      parallel = numprocs .gt. 1
      if(
     +  numprocs .gt. mnodes )call crash( 'GETSET', 'mnodes too small' )
c set all logical unit nos to nonsense values
      do i = 1, 6 + 2 * mgrids
        iunit( i ) = -1
      end do
c read run number from stdin
      call gtintg( 'Enter run no', irun )
      if( master )call opnfil( ni, 'dat', 'formatted', 'old', 'seq', i )
      call mpi_bcast( ni, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( i, 1, mpi_integer, 0, mpi_comm_world, ierr )
      if( i .ne. 0 )call crash( 'GETSET', '.dat file not found' )
c check run number
      call getline( ni, line )
      read( line( 11:80 ), * )i
      if( i .ne. irun )then
        if( master )then
          print *, 'Run number input was', irun
          print *, 'Number read from data file was', i
        end if
        call crash( 'GETSET', 'Incorrect data file for run' )
      end if
c open output file
      if( master
     +        )call opnfil( no, 'lis', 'formatted', 'old', 'append', i )
      call mpi_bcast( no, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( i, 1, mpi_integer, 0, mpi_comm_world, ierr )
      if( i .ne. 0 )then
        if( master )then
          print *, 'Old .lis file not found - opening a new file'
          call opnfil( no, 'lis', 'formatted', 'new', 'seq', i )
        end if
        call mpi_bcast( i, 1, mpi_integer, 0, mpi_comm_world, ierr )
        if( i .ne. 0 )call crash( 'GETSET', 'Error opening .lis file' )
      end if
c read and interpret field determination method
      call getline( ni, line )
      i = 11
      do while ( line( i:i ) .eq. ' ' )
        i = i + 1
      end do
      read( line( i:i+2 ), '( a )' )code
      hybrid = .false.
      twogrd = .false.
      if( ( code .eq. 'HYB' ) .or. ( code .eq. 'TWO' ) )then
        if( code .eq. 'HYB' )hybrid = .true.
        if( code .eq. 'TWO' )twogrd = .true.
        read( line( 21:80 ), * )ngrid
        call getline( ni, line )
        i = 11
        do while ( line( i:i ) .eq. ' ' )
          i = i + 1
        end do
        read( line( i:i+2 ), '( a )' )code
      else
        ngrid = 1
      end if
      do j = 1, ngrid
        call lowercase( code )
        if( master )print
     +        '( '' Field determination method selected: '', a3 )', code
        i = -1
        if( code .eq. 'non' )i = 0
        if( code .eq. 'p2d' )i = 1
        if( code .eq. 'c3d' )i = 2
        if( code .eq. 'sf2' )i = 3
        if( code .eq. 'p3a' )i = 4
        if( code .eq. 'p3d' )i = 5
        if( code .eq. 's3d' )i = 6
        if( code .eq. 'sf3' )i = 7
        if( code .eq. 'c2d' )i = 8
        if( code .eq. 'scf' )i = 9
        if( code .eq. 'dr3' )i = 10
        if( code .eq. 'bht' )i = 11
        if( i .lt. 0 )call crash( 'GETSET', 'Unrecognized method' )
        igrid( j ) = i
        if( j .lt. ngrid )then
          call getline( ni, line )
          i = 11
          do while ( line( i:i ) .eq. ' ' )
            i = i + 1
          end do
          read( line( i:i+2 ), '( a )' )code
        end if
      end do
      call codeset
c report selected option
      if( master )then
        if( ngrid .eq. 1 )then
          write( no,
     +        '( '' Field determination method selected: '', a3 )' )code
        else if( hybrid )then
          write( no, '( '' Hybrid methods selected: '', 3i2 )' )
     +                                      ( igrid( i ), i = 1, ngrid )
        else
          write( no, '( '' Multiple grids selected: '', 3i2 )' )
     +                                      ( igrid( i ), i = 1, ngrid )
        end if
      end if
      return
      end
