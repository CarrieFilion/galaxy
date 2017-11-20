      real*8 function povfun( t )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for Poveda's expression for the density in an R^{1/4} surface
c   density sphere
c
c calling argument
      real*8 t
c
c common block
c
      common / povedal / b, aj
      real*8 aj, b
c
      povfun = exp( -b * aj * t ) / sqrt( t**8 - 1. )
      return
      end
