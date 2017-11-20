      subroutine jshist( a, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to draw a histogram in the current plotting window of the
c   the n values in the array a
c   the number of histogram bins is suitably determined from n
c   the y-scale is set so that the total are under the histrogram is unity
c
c calling arguments
      integer n
      real a( n )
c
c local arrays
      integer mbin
      parameter ( mbin = 101 )
      integer bin( mbin )
c
c local variables
      integer i, j, nbin
      real dx, scale, x, xmax, xmin, y, ymax, ymin
c
c initialise bins
      x = n
      nbin = sqrt( x )
      nbin = min( nbin, mbin )
      do i = 1, nbin
        bin( i ) = 0
      end do
c set scale
      xmin = 1.e+10
      xmax = -xmin
      do i = 1, n
        xmax = max( xmax, a( i ) )
        xmin = min( xmin, a( i ) )
      end do
      dx = ( xmax - xmin ) / real( nbin )
      scale = 1. / dx
c bin up data
      do j = 1, n
        i = 1 + int( scale * ( a( j ) - xmin ) )
        i = min( i, nbin )
        bin( i ) = bin( i ) + 1
      end do
c set up frame
      call jspage
      j = 0
      do i = 1, nbin
        j = max( j, bin( i ) )
      end do
      scale = scale / real( n )
      ymax = 1.1 * scale * real( j )
      ymin = 0
      call jscale( xmin, xmax, ymin, ymax )
c normalise and plot
      x = xmin
      call jsmove( x, 0. )
      do i = 1, nbin
        y = scale * real( bin( i ) )
        call jsline( x, y )
        x = x + dx
        call jsline( x, y )
      end do
      call jsline( x, 0. )
c draw axes
      call jsaxis( 'x', 'a', 1 )
      call jsaxis( 'y', 'freq', 1 )
      return
      end
