      subroutine ncheck
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to perform a simple check that the number of particles contributing
c   to the active density, plus those both in the central hole and outside the
c   grid, is equal to the number in the simulation.
      use aarrays
      implicit none
c
c common blocks
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
      include 'inc/model.f'
c
      include 'mpif.h'
c
c local variables
      integer i, ig, j, jz, k, l
      logical fail, partial
      real*8 x, y
c
c check lists only if particles have differing masses or precision will be a problem
      if( uqmass .or. noslfg .or. dr3d .or. bht .or.
     +    ( nbod .gt. 1000000 ) )then
        call chklst
      else
c check impossible if a substep or for S3D when monopole terms are omitted
        if( lsfp )then
          partial = mod( istep, nstep( nzones ) ) .ne. 1
        else
          partial = mod( istep, nstep( nzones ) ) .ne. 0
        end if
        if( .not. ( partial .or. ( s3d .and. fixaxi ) ) )then
c sum mass on the various grids - couldn't I use grdmas here??
          jz = 0
          do ig = 1, ngrid
            call switch( ig )
            x = 0
            if( lsfp )then
              jz = jz + ncontrib
            else if( s3d )then
              l = 4 * s3ntm * ( nr( jgrid ) - 1 ) + 1
              x = x +
     +            s3dfld( 1, 1 ) + s3dfld( l, 1 ) * s3rad( nr( jgrid ) )
              jz = jz + nint( x )
            else if( lgrd )then
              do i = 1, mesh( jgrid )
                x = x + grdmss( i, 1 )
              end do
              if( p2d .or. p3d )then
c combine totals from different nodes for polar grids only
                call mpi_allreduce( x, y, 1,
     +                mpi_double_precision, mpi_sum, mpi_comm_world, i )
                x = y
              end if
              jz = jz + nint( x )
            end if
          end do
          call switch( 0 )
          if( master )then
c report if requested
    1       if( lprint )write( no, 204 )jz, mstep( 1 )
c check total number of particles
            k = 0
            do i = 1, 4
              k = k + ncen( i )
            end do
            if( lprint )then
              if( k .ne. 0 )write( no, 200 )k
              if( noff .ne. 0 )write( no, 201 )noff
              if( nshort .ne. 0 )write( no, 202 )nshort
            end if
            k = k + jz + noff
c allow for massless particles
            do j = 1, ncmp
              if( testp( j ) )k = k + nsp( j )
            end do
            if( s3d .and. abs( k - nbod ) .lt. 5 )k = nbod
            fail = k .ne. nbod
            if( fail )then
              if( lprint )write( no, 203 )x, mstep( 1 )
              lprint = .not. lprint
              if( lprint )go to 1
            end if
          end if
          call mpi_bcast( fail, 1, mpi_logical, 0, mpi_comm_world, j )
          if( fail )call crash( 'NCHECK',
     +                                 'Incorrect number of particles' )
        end if
      end if
      return
  200 format( i10, ' particles are in the central hole' )
  201 format( i10, ' particles are off the mesh' )
  202 format( i10, ' particles are on extra-short time steps' )
  203 format( f12.2, ' particles active at step', i7 )
  204 format( i10, ' particles active at step', i7 )
      end
