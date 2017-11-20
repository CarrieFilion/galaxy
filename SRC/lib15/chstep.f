      subroutine chstep( is, tfac )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine for changing the time step for a single particle.  It is needed
c   when a particle changes zones, or enters or leaves intensive care
c
c calling arguments
      integer is
      real tfac
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
c local variable
      integer i
c
      if( lfrv )then
c RVLF: advance the velocity half of the old step and back it up half of the new
        do i = 1, ndimen
          oldc( i + ndimen, is ) = tfac * ( oldc( i + ndimen, is ) +
     +                               .5 * acc( i, is ) * ( 1. - tfac ) )
          acc( i, is ) = acc( i, is ) * tfac * tfac
c step forward with the new time step
          newc( i + ndimen, is ) = oldc( i + ndimen, is ) + acc( i, is )
          newc( i, is ) = oldc( i, is ) + newc( i + ndimen, is )
        end do
      else
        call crash( 'CHSTEP', 'This option should never be needed' )
c APLF: back up the position half of the old step and advance it half of the new
        do i = 1, ndimen
          newc( i, is ) = newc( i, is ) +
     +                       .5 * ( tfac - 1. ) * newc( i + ndimen, is )
          newc( i + ndimen, is ) = tfac * newc( i + ndimen, is )
          acc( i, is ) = acc( i, is ) * tfac * tfac
        end do
      end if
      return
      end
