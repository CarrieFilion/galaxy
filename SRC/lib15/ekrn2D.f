      real function ekrn2D( xx, yy )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns simple kernel estimate of a 2-D function from discrete data
c
c calling arguments
      real xx, yy
c
c common blocks
c
      include 'inc/kernel.f'
c
c local variables
      integer i
      real h2, p, r2, x, y
      include 'inc/pi.f'
c
c check for sensible values only once
      if( krnuse .le. 0 )then
        if( hkrn .le. 0. )call crash( 'EKRN2D', 'Kernel width not set' )
        if( xper .and. ( xwidth .le. 0. ) )then
          print *, xwidth
          call crash( 'EKRN2D',
     +                       'periodic in x, but with nonsense period' )
        end if
        if( yper .and. ( ywidth .le. 0. ) )then
          print *, ywidth
          call crash( 'EKRN2D',
     +                       'periodic in y, but with nonsense period' )
        end if
        krnuse = 1
      end if
c initialize
      h2 = hkrn * hkrn
      ekrn2D = 0
c work over points
      do i = 1, nkrn
        p = akrn( 1, i )
        x = p - xx
        if( xper )then
          if( abs( x - xwidth ) .lt. hkrn )x = x - xwidth
          if( abs( x + xwidth ) .lt. hkrn )x = x + xwidth
        end if
        if( abs( x ) .lt. hkrn )then
          p = akrn( 2, i )
          y = p - yy
          if( yper )then
            if( abs( y - ywidth ) .lt. hkrn )y = y - ywidth
            if( abs( y + ywidth ) .lt. hkrn )y = y + ywidth
          end if
          if( abs( y ) .lt. hkrn )then
            r2 = x * x + y * y
c Epanechnikov kernel
            if( r2 .lt. h2 )ekrn2D = ekrn2D + 1. - r2 / h2
          end if
        end if
      end do
c normalize
      ekrn2D = 2. * ekrn2D / ( pi * h2 * real( nkrn ) )
      return
      end
