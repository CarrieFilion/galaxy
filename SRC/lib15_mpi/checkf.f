      subroutine checkf
c  Copyright (C) 2015, Jerry Sellwood
c
c generic routine to check solution for field of a few test masses on a grid
c   parallel version, but works in single processor mode too
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'mpif.h'
c
c local array
      integer stats( mpi_status_size )
c
c local variables
      integer i, ierr, inode, j, k, l, tag
      logical done, ok, okall, skip
c
      if( master )then
        done = .true.
        if( p2d )then
          call p2chkf( 1, nr( 1 ) )
        else if( p3d )then
          call p3chkf
        else if( c2d )then
          call c2chkf
        else if( p3a )then
          call pachkf
        else
          done = .false.
        end if
      end if
      call mpi_bcast( done, 1, mpi_logical, 0, mpi_comm_world, ierr )
      if( .not. done )call crash( 'CHECKF', 'Unrecognized grid' )
c
      if( parallel )then
        if( master )write( no, * )'Result on master node'
c compare solutions on other nodes
        l = mesh( jgrid )
        okall = .true.
        tag = 0
        k = 4
        if( twod )k = 3
        do j = 1, k
          skip = ( ( .not. phys ) .and. ( j .eq. k ) )
          if( .not. skip )then
            do inode = 1, numprocs - 1
              if( p2d .or. p3d )then
c copy results to master node
                if( myid .eq. inode )then
                  call mpi_send( grdfld( 1, j ), l, mpi_real, 0, tag,
     +                           mpi_comm_world, ierr )
                else if( myid .eq. 0 )then
                  call mpi_recv( grdmss( 1, 1 ), l, mpi_real, inode,
     +                           tag, mpi_comm_world, stats, ierr )
                end if
c compare
                if( master )then
                  ok = .true.
                  do i = 1, l
                    ok = ok .and. ( grdfld( i, j ) .eq. grdmss( i, 1 ) )
                  end do
                end if
              else
                call crash( 'CHECKF', 'grid type not programmed' )
              end if
              if( master )then
                if( ok )then
                  print *, 'OK for field type', j, ' from node', inode
                else
                  print *,
     +             'Check failed for field type', j, ' from node', inode
                  write( no, * )
     +             'Check failed for field type', j, ' from node', inode
                end if
                okall = okall .and. ok
              end if
            end do
          end if
        end do
c finish up
        if( master )then
          if( okall )then
            write( no, * )'Fields on all other nodes are identical'
          else
            write( no, * )'Fields differ on other nodes - watch out!'
          end if
        end if
      end if
      return
      end
