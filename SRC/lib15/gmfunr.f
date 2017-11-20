      real*8 function gmfunr( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for radius weighted mass integral, called initially from
c   supset
c
c calling argument
      real*8 r
c
c external
      real*8 gmfuni
c
      gmfunr = r * gmfuni( r )
      return
      end
