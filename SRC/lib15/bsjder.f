      real*8 function bsjder( m, x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns derivative wrt argument of J_m(x)
c
c calling arguments
      integer m
      real*8 x
c
c external
      real*8 bsjfun
c
c special arguments
      if( x .eq. 0.d0 )then
        if( m .eq. 1 )then
          bsjder = .5
        else
          bsjder = 0
        end if
      else if( m .eq. 0 )then
        bsjder = -bsjfun( 1, x )
      else
c other values - formula 9.1.27 of Abramowitz & Stegun
        bsjder = bsjfun( m - 1, x ) - m * bsjfun( m, x ) / x
      end if
      return
      end
