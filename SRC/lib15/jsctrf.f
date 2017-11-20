      subroutine jsctrf( funct, ncont, work, npx, npy, trans, pos )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to evaluate a given input function at an array of points
c   in the current plotting window and to draw ncont contours of them.
c   The contour levels are determined from the range of function values
c   found, with the logicals pos=.T. prescribing that only positive
c   values are to be contoured
c
c calling arguments
      external trans
      integer ncont, npx, npy
      logical pos
      real funct, work( npx, npy )
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
      real fmax, fmin, x, x1, x2, y, y1, y2
c
c set up scaling constants
      xscale = ( sx2 - sx1 ) / real( npx - 1 )
      yscale = ( sy2 - sy1 ) / real( npy - 1 )
      xoffs = sx1 - xscale
      yoffs = sy1 - yscale
c choose x and y
      do iy = 1, npy
        y1 = iy
        do ix = 1, npx
          x1 = ix
          call trans( x1, y1, x, y )
          work( ix, iy ) = funct( x, y )
        end do
      end do
c find minimum and maximum
      fmin = 1.e + 24
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
      if( fmax - fmin .lt. 1.e-7 )return
c mark these
      call jschsz( .1 )
      call trans( x2, y2, x, y )
      call jsymbl( x, y, 1, 5 )
      call jsbldt( 'max' )
      call jswrit( x, y )
      call trans( x1, y1, x, y )
      call jsymbl( x, y, 1, 5 )
      call jsbldt( 'min' )
      call jswrit( x, y )
      call jschsz( .3 )
c linear contour spacings
      if( pos )fmin = 0
      do ic = 1, ncont
        c( ic ) = fmin + ( real( ic - 1 ) + .5 ) * ( fmax - fmin )
     +                                                   / real( ncont )
      end do
c logarithmically spaced contours
c      fmax = log10( fmax )
c      do ic = 1, ncont
c        x = fmax - .4 * real( ic )
c        c( ic ) = 10.**x
c      end do
c contour function
      call jscont( work, npx, npy, c, ncont, trans )
      return
      end
