      real function ekrn1D( xx )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns simple kernel estimate of a 1-D function from discrete data
c
c calling argument
      real xx
c
c common blocks
c
      include 'inc/kernel.f'
c
c local variables
      integer i
      logical firstc
      real x
      save firstc
c
      data firstc / .true. /
c
c check for sensible values only once
      if( firstc )then
        if( hkrn .le. 0. )call crash( 'EKRN1D', 'Kernel width not set' )
        if( xper .and. ( xwidth .le. 0. ) )then
          print *, xwidth
          call crash( 'EKRN1D',
     +                       'periodic in x, but with nonsense period' )
        end if
        firstc = .false.
      end if
c work over points
      ekrn1D = 0
      do i = 1, nkrn
        x = akrn( 1, i ) - xx
        if( xper )then
          if( x - xwidth .gt. -hkrn )x = x - xwidth
          if( x + xwidth .lt. hkrn )x = x + xwidth
        end if
c Epanechnikov kernel
        if( abs( x ) .lt. hkrn )ekrn1D = ekrn1D + 1. - ( x / hkrn )**2
      end do
c normalize
      ekrn1D = .75 * ekrn1D / ( hkrn * real( nkrn ) )
      return
      end
