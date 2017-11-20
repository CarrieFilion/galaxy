      subroutine ptbscl( ready )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to rescale the acceleration of the perturber
c
c calling argument
      logical ready
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'mpif.h'
c
c external
      real tsfac
c
c local array
      real*8 abuf( 3, mzones )
c
c local variables
      integer i, j, jz, oldoff
      logical update
      real fac, fract
      save oldoff
c
c check whether a perturber is present
      if( pertbn )then
        if( ready )then
c sum results from all processors if required
          if( parallel )then
            call mpi_allreduce( nptrb, j, 1, mpi_integer, mpi_sum,
     +                                               mpi_comm_world, i )
            nptrb = j
            call mpi_allreduce( pzacc, abuf, 3 * mzone,
     +                mpi_double_precision, mpi_sum, mpi_comm_world, i )
            do jz = 1, mzone
              do i = 1, 3
                pzacc( i, jz ) = abuf( i, jz )
              end do
            end do
          end if
c check particle counter of accelerations from perturber
          if( mod( istep, nstep( nzones ) ) .eq. 0 )then
            nptrb = nptrb + oldoff
            if( nptrb .ne. nbod )then
              if( master )print '( i10, '' particles at step'', i8 )',
     +                                                      nptrb, istep
              call crash( 'PTBSCL', 'Wrong nptrb' )
            end if
          end if
c save newly computed values and discard obsolete ones only once per step
          update = mzone .eq. 1
          if( .not. update )update =
     +                             mstep( mzone ) .eq. lstep( 2, mzone )
          if( update )then
            do jz = 1, mzone
              do i = 1, 3
                accpz( i, 1, jz ) = accpz( i, 2, jz )
                accpz( i, 2, jz ) = pzacc( i, jz )
              end do
            end do
          else
            call crash( 'PTBSCL', 'Acceleration update error' )
          end if
c update current acceleration
          do i = 1, 3
            accptb( i ) = accpz( i, 2, 1 )
          end do
c combine accelerations from different zones
          if( nzones .gt. 1 )then
            do jz = 2, nzones
c interpolate linearly between times
              j = lstep( 2, jz ) - istep
c            fac = 1
c            if( .not. quapp )
              fac = 1. / tsfac( jz )**2
              if( j .ne. nstep( jz ) )then
                fract = 1. - real( j ) / real( nstep( jz ) )
c                print 200, fract, jz, lstep( 2, jz )
c  200 format( ' Adding fraction', f5.2, ' of accels from zone no', i2,
c     +        ' step', i7 )
                do i = 1, 3
                  accptb( i ) = accptb( i ) +
     +                                   fract * accpz( i, 2, jz ) * fac
                end do
              else
                fract = 0.
              end if
              if( j .ne. 0 )then
                fract = 1. - fract
c                print 200, fract, jz, lstep( 1, jz )
                do i = 1, 3
                  accptb( i ) = accptb( i ) +
     +                                   fract * accpz( i, 1, jz ) * fac
                end do
              end if
            end do
          end if
          if( exsphr .or. exgnrc )then
c re-scale acceleration of perturber by ratio of particle to perturber masses
            fac = pmass / ( lscale**3 * ts**2 * mptbr )
            do i = 1, 3
              accptb( i ) = fac * accptb( i )
            end do
          else if( extbar )then
c convert torque to external units
c          if( .not. quapp )
            accptb( 1 ) = accptb( 1 ) / ( lscale * ts )**2
c scale torque on perturber by particle mass and switch sign for back reaction
            fac = pmass / ( lscale**3 * ts**2 )
            accptb( 1 ) = -fac * accptb( 1 )
c acceleration of the bar in response to applied torque
c    MoI for a homogeneous ellipsoid is M(a^2+b^2)/5
            accptb( 1 ) = 5. * accptb( 1 ) /
     +                                 ( mptbr * ( abar**2 + bbar**2 ) )
c no acceleration of imposed spiral
          else if( exspir )then
            do i = 1, 3
              accptb( i ) = 0
            end do
          else
            call crash( 'PTBSCL', 'Unknown perturber' )
          end if
        else
c initialize acceleration components of perturber
          sumptb = .true.
          do jz = 1, mzone
            do i = 1, 3
              pzacc( i, jz ) = 0
            end do
          end do
          if( phys )peptrb = 0
          nptrb = 0
          oldoff = noff
        end if
      end if
      return
      end
