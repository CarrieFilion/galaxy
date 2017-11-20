      subroutine p3manlp( jrs, krs )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to perform the double Fourier analysis of the mass array for the
c   3-D polar grid.  The Fourier coefficients overwrite the input array
c   and are double the size
c
c This version uses single precision FFTPAK routines
      use aarrays
      implicit none
c
c calling arguments
      integer jrs, krs
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
      integer ltrig
      real, allocatable :: trig(:)
c
      integer lwork
      real, allocatable :: w(:)
c
c local variables
      integer ia, ierr, im, ir, ispac, iz, j, m, mtrig, nm, nrg, nrs
      real strig
c
      nrg = nr( jgrid )
c determine size of transformed active mass region
      nrs = krs - jrs + 1
      nm = 0
      do im = 1, ng
        if( .not. lg( im ) )then
          if( im .eq. 1 .or. im .eq. mmax )then
            nm = nm + 1
          else
            nm = nm + 2
          end if
        end if
      end do
      ispac = nrs * ngz * nm
c allocate space
      lwork = 4 * ispac
      allocate ( w( lwork ) )
      ltrig = 2 * max( na, 2 * ngz ) + 15
      allocate ( trig( ltrig ) )
c initialize if the vector length is new
      if( na .ne. mtrig )then
c check space
        mtrig = na
        j = 2 * mtrig + 15
        if( j .gt. ltrig )call space( ltrig, j, '/ trig /', 'P2MANLP' )
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        strig = 1. / real( mtrig )
      end if
c Fourier analysis in azimuth of each grid ring and plane in turn
      j = 1
c      ic = 0
      do ir = jrs, krs
        do iz = 1, ngz
c no parallelization here
c          if( mod( ic, numprocs ) .eq. myid )then
            call sfftf1( na, grdmss( j, 0 ),
     +                        trig, trig( na + 1 ), trig( 2 * na + 1 ) )
c          else
c            do ia = j, j + na - 1
c              msh( ia ) = 0
c            end do
c          end if
c          ic = ic + 1
          j = j + na
        end do
c        j = j + ( nrg - nrs ) * ngz * na
      end do
c rearrange and double vertical size of mass array
      m = 0
      do ir = 1, nrs
        do ia = 1, na
c copy only selected sectoral harmonics: im = m + 1
          im = ia / 2 + 1
          if( .not. lg( im ) )then
            j = ( ir - 1 ) * ngz * na + ia
            do iz = 1, ngz
              m = m + 1
              w( m ) = grdmss( j, 0 )
              w( m + ngz ) = 0
              j = j + na
            end do
            m = m + ngz
          end if
        end do
      end do
c this operation combines masses from all the grids
      call mpi_allreduce( w, grdmss( 1, 0 ), m, mpi_real,
     +                    mpi_sum, mpi_comm_world, ierr )
c initialize if the vector length is new
      if( 2 * ngz .ne. mtrig )then
        mtrig = 2 * ngz
c check space
        j = 2 * mtrig + 15
        if( j .gt. ltrig )call space( ltrig, j, '/ trig /', 'PAMANLP' )
        call sffti1( mtrig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
        strig = 1. / real( mtrig )
      end if
c Fourier analysis of each vector in turn - not parallelized
      j = 1
      do m = 1, nrs * nm
c        if( mod( m - 1, numprocs ) .eq. myid )then
          call sfftf1( mtrig, grdmss( j, 0 ),
     +                  trig, trig( mtrig + 1 ), trig( 2 * mtrig + 1 ) )
c        else
c          do iz = j, j + mtrig - 1
c            grdmss( iz, 0 ) = 0
c          end do
c        end if
        j = j + mtrig
      end do
c      call mpi_allreduce( grdmss( 1, 0 ), w, 2 * ispac,
c     +                    mpi_real, mpi_sum, mpi_comm_world, ierr )
c copy back results
c      do j = 1, 2 * ispac
c        grdmss( j, 0 ) = w( j )
c      end do
      return
      end
