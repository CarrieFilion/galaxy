      real*8 function frext( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the central attraction of the exteral mass added into
c   an active one to compress its DF
c
c calling argument
      real*8 r
c
c external
      real*8 gmext
c
c local variable
      real*8 r1
c
      r1 = max( abs( r ), 1.d-8 )
      frext = -gmext( r1 ) / r1**2
      return
      end
