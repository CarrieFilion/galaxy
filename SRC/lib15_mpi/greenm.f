      subroutine greenm( direct )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c driver routine to create the data tables needed for the field determination
c
c Grid methods generally require a Greens function.  The routine passes on
c   the calling argument, which indicates whether or not a direct version is
c   required in addition to the default transformed one.
c
c SFP methods need tables of function values and their derivatives.
c
c No tables are required for the PM+SH method.
c
c calling argument
      logical direct
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'mpif.h'
c
c local variables
      integer i, j( 4 ), k
c
      data j / 4 * 0 /
c
      do i = 1, ngrid
        call switch( i )
        if( ncode .gt. 11 )call crash( 'GREENM', 'Unrecognized method' )
c grid methods only
        if( c3d )then
          call crash( 'GREENM', 'Parallel version not available' )
        else if( p2d .or. p3d .or. c2d .or. p3a )then
c open and check the file(s) on the master node only
          if( master )then
            call opnfil( ngt, 'grt', 'unformatted', 'old', 'seq', j )
            if( j( 1 ) .eq. 0 )then
              call grnhed( .true., .false., j( 2 ) )
              if( ( j( 2 ) .eq. 0 ) .and. direct )then
                call opnfil( ngd,
     +                      'grd', 'unformatted', 'old', 'seq', j( 3 ) )
c verify existing direct file
                if( j( 3 ) .eq. 0 )then
                  call grnhed( .true., .true., j( 4 ) )
                else
                  j( 4 ) = 0
                end if
              end if
            end if
          end if
c broadcast results
          call mpi_bcast( j, 4, mpi_integer, 0, mpi_comm_world, k )
          if( j( 2 ) .ne. 0 )call crash( 'GREENM',
     +                             'Wrong transformed Green fun file' )
          if( direct )then
            if( j( 3 ) .ne. 0 )call crash( 'GREENM',
     +                               'Direct Green fn file not found' )
            if( j( 4 ) .ne. 0 )call crash( 'GREENM',
     +                                  'Wrong direct Green fun file' )
          end if
c make the file(s) on the master node only
          if( j( 1 ) .ne. 0 )then
            if( master )then
              if( p2d )call p2grnm( direct )
              if( c3d )call c3grnm
              if( p3d )call p3grnm( direct )
              if( p3a )call pagrnm( direct )
              if( c2d )call c2grnm( direct )
              print *, 'Green fun file created for grid', jgrid
            end if
c force slaves to wait until master has finished
            call mpi_barrier( mpi_comm_world, k )
          end if
          if( .not. master )then
c open the .grt file only on the slaves - it is already opened on the master
            call opnfil( ngt, 'grt', 'unformatted', 'old', 'seq', k )
          end if
        end if
      end do
      call switch( 0 )
      return
      end
