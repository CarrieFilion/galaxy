      subroutine massum
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c Routine to combine mass arrays from different nodes, when needed
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'mpif.h'
c
c local allocatable arrays
      real, allocatable :: w(:,:)
      real, allocatable :: x(:,:)
c
c local array
      integer irad( mzones )
c
c local variables
      integer i, ig, j
      integer nzone
c
c do this only if parallel is on
      if( parallel )then
c work over active grids
        do ig = 1, ngrid
          call switch( ig )
          if( dr3d )then
c copy and combine new positions
            allocate ( w( 3, ndrct ) )
            allocate ( x( 3, ndrct ) )
            do i = 1, ndrct
              do j = 1, 3
                w( j, i ) = drpt( j + 4, i )
              end do
            end do
            call mpi_allreduce( w, x, 3 * ndrct,
     +                          mpi_real, mpi_sum, mpi_comm_world, i )
c update old positions
            do i = 1, ndrct
              do j = 1, 3
                drpt( j + 1, i ) = x( j, i )
              end do
            end do
            deallocate ( w )
            deallocate ( x )
          else if( .not. noslfg )then
c grid or field methods
            if( p2d .or. p3d )then
c need only reset radial range markers - masses are combined later
              call mpi_allreduce(
     +     jrad, irad, nzones, mpi_integer, mpi_min, mpi_comm_world, i )
              do i = 1, nzones
                jrad( i ) = irad( i )
              end do
              call mpi_allreduce(
     +     krad, irad, nzones, mpi_integer, mpi_max, mpi_comm_world, i )
              do i = 1, nzones
                krad( i ) = irad( i )
              end do
            else
c use temporary work space - double size for real*8 data
              j = mesh( jgrid )
              if( lsfp )then
                allocate ( w( j, 2 ) )
              else
                allocate ( w( j, 1 ) )
              end if
              do nzone = 1, mzone
                if( lsfp )then
                  call mpi_allreduce( sfpmss( 1, nzone ), w, j,
     +                mpi_double_precision, mpi_sum, mpi_comm_world, i )
                  call blkcpy2( w, sfpmss( 1, nzone ), j )
                else if( s3d )then
                  call mpi_allreduce( s3dmss( 1, nzone, jgrid ), w, j,
     +                            mpi_real, mpi_sum, mpi_comm_world, i )
                  call blkcpy( w, s3dmss( 1, nzone, jgrid ), j )
                  if( hybrid .and. ( jgrid .eq. 2 ) )then
                    call mpi_allreduce( s3dmss( 1, nzone, 1 ), w, j,
     +                            mpi_real, mpi_sum, mpi_comm_world, i )
                    call blkcpy( w, s3dmss( 1, nzone, 1 ), j )
                  end if
                else
                  call mpi_allreduce( grdmss( 1, nzone ), w, j,
     +                            mpi_real, mpi_sum, mpi_comm_world, i )
                  call blkcpy( w, grdmss( 1, nzone ), j )
                end if
              end do
              deallocate ( w )
            end if
          end if
        end do
      end if
      call switch( 0 )
      return
      end
