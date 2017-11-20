      real*8 function celli1( arg )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      real*8 arg
c
c external
      real*8 celint
c
c local variables
      real*8 aa, bb, kc, pp
      parameter ( aa = 1, bb = 1, pp = 1 )
c
      kc = sqrt( 1. - arg * arg )
      kc = max( kc, 1.d-20 )
      celli1 = celint( kc, pp, aa, bb )
      return
      end
