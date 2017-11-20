      real*8 function celli2( arg )
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
      parameter ( aa = 1, pp = 1 )
c
      bb = 1. - arg * arg
      kc = sqrt( bb )
      celli2 = celint( kc, pp, aa, bb )
      return
      end
