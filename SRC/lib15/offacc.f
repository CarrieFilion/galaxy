      subroutine offacc( jst )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c adds the forces from particles off the grid - those that leave the grid
c   at this step should not be included
c
c calling arguments
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local array
      real xd( 3 ), xoff( 3 )
c
c local variables
      integer i, is, j, js
      logical notfound
      real a, d2, d3
c
c skip if nothing to do
      if( noff .gt. 0 )then
c assume no particles are off the grid during set up
        if( istep .lt. 0 )call crash( 'OFFACC',
     +                      'Some particles are off the grid at start' )
        if( parallel )call crash( 'OFFACC', ' Parallel version needed' )
c if particles in this list are on the grid
        if( ilist .lt. nlists )then
c work over particles off grid
          do js = 1, noff
            j = loco( js )
            if( j .ge. 0 )then
              do i = 1, ndimen
                xoff( i ) = ptcls( j + i )
              end do
c work over particles in this group
              do is = 1, jst
c distance between particles
                d2 = softl2
                do i = 1, ndimen
                  xd( i ) = oldc( i, is ) - xoff( i )
                  d2 = d2 + xd( i )**2
                end do
                d3 = d2**1.5
c acceleration components from this off-grid particle
                do i = 1, ndimen
                  a = -pmass * xd( i ) / d3
                  acc( i, is ) = acc( i, is ) + a
                  accoff( i, js ) = accoff( i, js ) - a
                end do
c sum PE if required
                if( phys )then
                  gpot( is ) = gpot( is ) - pmass / sqrt( d2 )
                  potoff( js ) = potoff( js ) - pmass / sqrt( d2 )
                end if
              end do
            end if
          end do
        else
c particles off the grid - accelerations now ready
          do is = 1, jst
            notfound = .true.
            do js = 1, noff
              if( loc( is ) .eq. loco( js ) )then
                notfound = .false.
                do i = 1, ndimen
                  acc( i, is ) = acc( i, is ) + accoff( i, js )
                end do
                gpot( is ) = gpot( is ) + potoff( js )
              end if
            end do
            if( notfound )then
              print *, noff, jst, is, ilist
              print *, ( loc( i ), i = 1, jst )
              print *, ( loco( js ), js = 1, noff )
              call crash( 'OFFACC', 'Particle not matched' )
            end if
          end do
c          print *, 'All particles matched in OFFACC', noff, jst
c          print *, ( loc( i ), i = 1, jst )
        end if
      end if
      return
      end
