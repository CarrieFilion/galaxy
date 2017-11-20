      subroutine rskern( before )
c routine to rescale the data values before contouring a 2D set of points
c   in order to make the kernel round.  The value of hkrn is assumed to be
c   that for the x-coordinate values
c a call with argument .false., which should not be made before the contouring
c   is completed, restores the plotting window to the scaling that was reset
c   in the first call
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      logical before
c
c common blocks
c
      include 'inc/jscmmn.f'
c
      include 'inc/kernel.f'
c
c local variables
      integer i
      real ar, ax1, ax2, ay1, ay2, xs, ys
      save ar, ax1, ax2, ay1, ay2
c
      if( before )then
c determine aspect ratio of the current plotting window
        ar = ( fy2 - fy1 ) / ( fx2 - fx1 )
c scaling factors - new x-range will be 0 -> 1 and y-range will be 0 -> ar
        xs = 1. / ( sx2 - sx1 )
        ys = ar / ( sy2 - sy1 )
        do i = 1, nkrn
          akrn( 1, i ) = xs * ( akrn( 1, i ) - sx1 )
          akrn( 2, i ) = ys * ( akrn( 2, i ) - sy1 )
        end do
c rescale kernel width
        hkrn = xs * hkrn
c preseve old scaling, then rescale the window
        ax1 = sx1
        ax2 = sx2
        ay1 = sy1
        ay2 = sy2
        call jscale( 0., 1., 0., ar )
      else
c restore orignal scaling
        hkrn = hkrn * ( ax2 - ax1 )
        call jscale( ax1, ax2, ay1, ay2 )
      end if
      return
      end
