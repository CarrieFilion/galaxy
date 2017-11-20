      subroutine p2manlp( ip, jr, kr )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to perform the Fourier analysis of the mass array for the 2-D
c   polar grid.  The Fourier coefficients overwrite the input array
c the calling paramaters are the base address of the array containing
c   the input data and the inner and outer radii of the partial grid
c   that is active (could also be the boundaries)
      use aarrays
      implicit none
c
c calling arguments
      integer ip, jr, kr
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
      real, allocatable :: w(:)
c
c local variables
      integer im, ir, j, n, nm, nrngs, pmesh
c
      ltrig = 2 * na + 15
      allocate ( trig( ltrig ) )
c initialize trig array
      call sffti1( na, trig( na + 1 ), trig( 2 * na + 1 ) )
c determine base address and size of active mass region
      j = ( jr - 1 ) * na + 1
c Fourier analysis of each grid ring in turn
      do ir = jr, kr
        call sfftf1( na, grdmss( j, ip ),
     +                        trig, trig( na + 1 ), trig( 2 * na + 1 ) )
        j = j + na
      end do
      deallocate ( trig )
c
      if( .not. parallel )call crash( 'P2FNDFP',
     +    'Parallel version called for a single processor calculation' )
c count active sectoral harmonics
      nm = 0
      do im = 1, mmax
        if( .not. lg( im ) )then
          if( im .eq. 1 .or. im .eq. mmax )then
            nm = nm + 1
          else
            nm = nm + 2
          end if
        end if
      end do
      nrngs = kr - jr + 1
      pmesh = nrngs * nm
      allocate ( w( pmesh ) )
c assemble only needed results
      n = ( jr - 1 ) * na
      do ir = jr, kr
        j = ( ir - jr ) * nm
        do im = 1, mmax
          n = n + 1
          if( .not. lg( im ) )then
            j = j + 1
            w( j ) = grdmss( n, ip )
          end if
          if( im .gt. 1 .and. im .lt. mmax )then
            n = n + 1
            if( .not. lg( im ) )then
              j = j + 1
              w( j ) = grdmss( n, ip )
            end if
          end if
        end do
      end do
c combine results from different processors
      n = ( jr - 1 ) * na + 1
      call mpi_allreduce(
     + w, grdmss( n, ip ), pmesh, mpi_real, mpi_sum, mpi_comm_world, j )
      call blkcpy( grdmss( n, ip ), w, pmesh )
c restore combined results
      n = ( jr - 1 ) * na
      j = 0
      do ir = jr, kr
        do im = 1, mmax
          n = n + 1
          if( lg( im ) )then
            grdmss( n, ip ) = 0
          else
            j = j + 1
            grdmss( n, ip ) = w( j )
          end if
          if( im .gt. 1 .and. im .lt. mmax )then
            n = n + 1
            if( lg( im ) )then
              grdmss( n, ip ) = 0
            else
              j = j + 1
              grdmss( n, ip ) = w( j )
            end if
          end if
        end do
      end do
      deallocate ( w )
      return
      end
