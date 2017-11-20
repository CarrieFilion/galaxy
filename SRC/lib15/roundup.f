      real function roundup( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns a value that is rounded up from the input value to some
c   (hopefully) intelligent value
c
c calling argument
      real x
c
c local array
      real tab( 10 )
c
c local variables
      integer i, j
      real t
c
      data tab / 1., 2., 3., 4., 5., 6., 8., 10., 10., 10. /
c
      roundup = x
      if( x .ne. 0. )then
        t = abs( x )
        i = log10( t )
        if( t .lt. 1. )i = i - 1
        j = t / 10.**i + 1.
        j = min( j, 10 )
        roundup = tab( j ) * 10.**i
        roundup = sign( roundup, x )
      end if
      return
      end
