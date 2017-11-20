      subroutine colplt( work, mx, my )
c  Copyright (C) 2014, Jerry Sellwood
c
c routine to make a color image of the values in array work - usually
c   relative overdensity
      implicit none
c
c calling arguments
      integer mx, my
      real work( mx, my )
c
c common block
c
      include 'inc/jscmmn.f'
c
c external
      logical gtlogl
c
c logan arrays
      real tr( 6 ), bl( 9 ), br( 9 ), bg( 9 ), bb( 9 )
c
c local variables
      integer i, icall, j
      logical logs
      real bright, contra, fmax, fmin, fmin2, fmin3
      parameter ( contra = 1, bright = 0.5 )
      save icall, logs, tr
c
      data icall / 0 /
c modified rainbow from white to black
      data bl / 0.0, 0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95, 1.0 /
      data br / 1.0, 0.5,  0.0,  0.0,  0.0,  1.0,  1.0,  0.5,  0.0 /
      data bg / 1.0, 0.0,  0.0,  1.0,  0.8,  0.8,  0.0,  0.0,  0.0 /
      data bb / 1.0, 0.5,  1.0,  1.0,  0.1,  0.0,  0.0,  0.0,  0.0 /
      data tr / 0., 0., 0., 0., 0., 0. /
c
      if( icall .eq. 0 )then
c initialize plot
        logs = gtlogl( 'logarithmic scale' )
        call pgctab( bl, br, bg, bb, 9, contra, bright )
c set up spatial scaling
        tr( 2 ) = ( sx2 - sx1 ) / real( mx - 1 )
        tr( 6 ) = ( sy2 - sy1 ) / real( my - 1 )
        tr( 1 ) = sx1 - tr( 2 )
        tr( 4 ) = sy1 - tr( 6 )
        icall = 1
      end if
c set density scale
      if( logs )then
        fmin = -5
        fmax = .2

        fmin = -1.1
        fmin2 = -2.4
        fmin3 = -3.7
        do j = 1, my
          do i = 1, mx
            if( work( i, j ) .gt. -50.
     +           .and. work( i, j ) .lt. fmin3  )then
              work( i, j ) = work( i, j ) + 3.9
            else if( work( i, j ) .lt. fmin2 )then
              work( i, j ) = work( i, j ) + 2.6
            else if( work( i, j ) .lt. fmin )then
              work( i, j ) = work( i, j ) + 1.3
            end if
          end do
        end do

      else
        fmax = .5
        fmin = -fmax
      end if
c output image
      call pgimag( work, mx, my, 1, mx, 1, my, fmin, fmax, tr )
      if( ipic .eq. 1 )call pgwedg( 'BI', 0.0, 3.0, fmin, fmax, ' ' )
      return
      end
