      real*8 function qtoomi( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes Q from a double integration over distribution function
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/moment.f'
c
c externals
      real*8 akappa, velint
c
c local variables
      real*8 gsig, sigu
c
c set power factors for velocity weighting
      iu = 2
      iv = 0
      sigu = velint( r )
      qtoomi = 0
      if( sigu .gt. 0.d0 )then
c normalise
        iu = 0
        iv = 0
        gsig = velint( r )
        if( gsig .gt. 0.d0 )then
          sigu = sqrt( sigu / gsig )
        else
          sigu = 1.d3
        end if
c compute q
        qtoomi = sigu * akappa( r ) / ( 3.358 * gsig )
        qtoomi = min( qtoomi, 10.d0 )
      end if
      return
      end
