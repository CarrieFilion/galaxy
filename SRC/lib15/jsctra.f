      subroutine jsctra( ncont, work, npx, npy, trans, pos, ln )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw ncont contours of the values in the input array
c   work.  The contour levels are determined from the range of values
c   in the input array, with the logicals pos=.T. prescribes that only
c   positive values to be contoured and ln=.T. prescribes that the
c   levels are spaced logarithmically, or linearly otherwise
c
c calling arguments
      external trans
      integer ncont, npx, npy
      logical ln, pos
      real work( npx, npy )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      include 'inc/jscmmn.f'
c
c local array
      real c( 50 )
c
c local variables
      integer ic, ix, iy
      real fmax, fmin, size, x, x1, x2, y, y1, y2
c
c set up scaling constants
      xscale = ( sx2 - sx1 ) / real( npx - 1 )
      yscale = ( sy2 - sy1 ) / real( npy - 1 )
      xoffs = sx1 - xscale
      yoffs = sy1 - yscale
c find minimum and maximum
      fmin = 1.e+24
      fmax = -fmin
      do iy = 1, npy
        do ix = 1, npx
          if( work( ix, iy ) .lt. fmin )then
            fmin = work( ix, iy )
            x1 = ix
            y1 = iy
          end if
          if( work( ix, iy ) .gt. fmax )then
            fmax = work( ix, iy )
            x2 = ix
            y2 = iy
          end if
        end do
      end do
c mark these
      size = height
      call jschsz( .1 )
      call trans( x2, y2, x, y )
      call jsymbl( x, y, 1, 5 )
      call jsbldt( 'max' )
      call jswrit( x, y )
      call trans( x1, y1, x, y )
      call jsymbl( x, y, 1, 5 )
      call jsbldt( 'min' )
      call jswrit( x, y )
      call jschsz( size )
c contour spacings
      if( ln )then
        if( fmax .le. 0. )then
          print *, 'JSCTRA called for logarithmic contours'
          print *, 'No positive values in array - fmax =', fmax
          return
        end if
c        fmax = max( fmax, 1.e-12 )
        fmax = log10( fmax )
c        fmin = max( fmin, 1.e-12 )
c        fmin = log( fmin )
        fmin = fmax - .4 * real( ncont )
      end if
      if( pos )fmin = 0
      do ic = 1, ncont
        c( ic ) = fmin + ( real( ic - 1 ) + .5 ) * ( fmax - fmin )
     +                                                   / real( ncont )
        if( ln )c( ic ) = 10.**c( ic )
      end do
c contour function
      call jscont( work, npx, npy, c, ncont, trans )
      return
      end
