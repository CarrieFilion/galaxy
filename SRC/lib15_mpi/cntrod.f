      subroutine cntrod
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to work through all the particles to obtain an improved estimate
c   of the position of the particle centroid on each grid.  Uses repeated
c   Newton-Raphson iterations until the change is acceptably small.
c
c Uses McGlynn's definition of the centroid (xp,yp,zp) that minimizes
c
c          omega_k = Sum [(x-xp)^2 + (y-yp)^2 + (z-zp^2)]^k
c
c   k=1 is equivalent to center of mass, which gives high weight to distant
c   particles.  I prefer k=1/2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      real w( 1 )
      include 'inc/buffer.f'
      equivalence ( w( 1 ), oldc( 1, 1 ) )
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
      integer i, is, iter, j, jst, k
      logical mp
      real a, b, c( 3 ), fx( 3, mgrids ), fxp( 3, mgrids )
c
c nothing to do if option is turned off or not a centering step
      if( ltroid )then
        if( kci .ne. 1 )call crash( 'CNTROD', 'Option not programmed' )
c iterate until shift is less than 0.1 percent of lscale
        a = lscale
        iter = 0
        do while ( ( a .gt. .001 * lscale ) .and. ( iter .lt. 4 ) )
          iter = iter + 1
c initialize
          do j = 1, ngrid
            do i = 1, 3
              fx( i, j ) = 0
              fxp( i, j ) = 0
            end do
          end do
c work through all lists except that of particles off the grid
          do ilist = 1, nlists - 1
c set izone & jlist
            call interpret
c work through groups of particles
            inext = islist( 1, ilist, myid + 1 )
            jgrid = jlist
            do while ( inext .ge. 0 )
              call gather( jst )
c work through this group of particles
              do is = 1, jst
                mp = cmpmas( iflag( is ) ) .gt. 0.
c skip massless particles
                if( mp )then
                  a = 0
                  do j = 1, 3
                    c( j ) = oldc( j, is ) - xcen( j, jgrid )
                    a = a + c( j )**2
                  end do
                  a = sqrt( a )
                  do j = 1, 3
                    fx( j, jgrid ) = fx( j, jgrid ) - c( j ) / a
                    b = 0
                    do k = 1, 3
                      if( j .ne. k )b = b + c( k )**2
                    end do
                    fxp( j, jgrid ) = fxp( j, jgrid ) + b / a**3
                  end do
                end if
              end do
            end do
          end do
c combine data from different processors
          if( parallel )then
c copy fx to a temporary buffer and sum
            k = 0
            do j = 1, ngrid
              do i = 1, 3
                k = k + 1
                w( k ) = fx( i, j )
              end do
            end do
            call mpi_allreduce( w, fx, k, mpi_real, mpi_sum,
     +                          mpi_comm_world, i )
c copy fxp to a temporary buffer and sum
            k = 0
            do j = 1, ngrid
              do i = 1, 3
                k = k + 1
                w( k ) = fxp( i, j )
              end do
            end do
            call mpi_allreduce( w, fxp, k, mpi_real, mpi_sum,
     +                          mpi_comm_world, i )
          end if
c combine results to keep centers coincident for hybrid methods
          mp = .false.
          if( ngrid .gt. 1 )then
            if( hybrid .and. ( ngrid .ne. 2 ) )call crash( 'CNTROD',
     +                                               'Whats going on?' )
            mp = hybrid .or. heavies
          end if
          if( mp )then
            do i = 1, 3
              fx( i, 1 ) = fx( i, 1 ) + fx( i, 2 )
              fx( i, 2 ) = fx( i, 1 )
              fxp( i, 1 ) = fxp( i, 1 ) + fxp( i, 2 )
              fxp( i, 2 ) = fxp( i, 1 )
            end do
          end if
c revise xcen
          a = 0
          do j = 1, ngrid
            b = 0
            do i = 1, 3
              b = b + ( fx( i, j ) / fxp( i, j ) )**2
              xcen( i, j ) = xcen( i, j ) - fx( i, j ) / fxp( i, j )
            end do
            b = sqrt( b )
            if( lprint )write( no, * )'center shift for grid', j, b
            a = max( a, b )
          end do
        end do
      else
        call crash( 'CNTROD', 'unknown rule' )
      end if
      return
      end
