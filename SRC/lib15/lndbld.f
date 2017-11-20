      real*8 function lndbld( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the value of l*kappa + m*Omega at the radius r
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/orbanl.f'
c
c externals
      real*8 akappa, omegac
c
      lndbld = real( morb ) * omegac( r )
      if( lorb .ne. 0 )lndbld = lndbld + real( lorb ) * akappa( r )
      return
      end
