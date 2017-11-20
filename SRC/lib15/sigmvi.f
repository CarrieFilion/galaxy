      real*8 function sigmvi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes sigma_v from a double integral over distribution fn
c
c calling argument
      real*8 r
c
c externals
      real*8 vmeani, v2meani
c
      sigmvi = v2meani( r ) - vmeani( r )**2
      sigmvi = max( sigmvi, 0.d0 )
      sigmvi = sqrt( sigmvi )
      return
      end
