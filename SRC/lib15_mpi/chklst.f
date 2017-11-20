      subroutine chklst
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to count particles in lists and verify totals, a substitute for
c   ncheck when individual particle masses mean that the number of active
c   particles cannot be determined from the mass assigned to the grid
c
c common blocks
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
      integer i, ia, ierr, j, k, m, n
      logical fail
      real a
      equivalence ( a, ia )
c
    1 n = 0
c work through lists
      do i = 1, nlists
        j = islist( 1, i, myid + 1 )
        k = 0
c count particles
        do while ( j .ge. 0 )
          k = k + 1
          a = ptcls( j + nwpp )
          j = ia
        end do
        if( parallel )then
          call mpi_allreduce( k, m, 1, mpi_integer,
     +                        mpi_sum, mpi_comm_world, ierr )
          k = m
        end if
        if( master .and.
     +      lprint )write( no,'( i8, '' particles in list'', i4 )' )k, i
        n = n + k
      end do
      if( master .and.
     +    lprint )write( no, '( i8, '' particles total'' )' )n
c checks
      call mpi_bcast( noff, 1, mpi_integer, 0, mpi_comm_world, i )
      fail = k .ne. noff
      if( fail )then
        if( master )print *, k, ' particles in list', nlists,
     +                                    ', whereas', noff, ' expected'
        call crash( 'CHKLST', 'Wrong number off grid' )
      end if
      call mpi_bcast( n, 1, mpi_integer, 0, mpi_comm_world, i )
      fail = n .ne. nbod
      if( fail )then
        if( master )print *, n, ' particles in all lists, whereas',
     +                                                nbod,  ' expected'
        lprint = .not. lprint
        if( lprint )go to 1
        call crash( 'CHKLST', 'Wrong total' )
      end if
      return
      end
