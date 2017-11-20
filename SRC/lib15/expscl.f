      subroutine expscl
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c multiplies the complex amplitude and phase data from a simulation
c   by an exponential factor that decreases with time in order to
c   make the fitted data more nearly equal in amplitude
c part of the mode fitting software
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
c local variables
      integer ip, it, n
      real scale
c
      if( apodc .eq. 0. )return
      n = 0
      do it = jt, kt
        scale = exp( apodc * ( tme( kt ) - tme( it ) ) )
        do ip = jp, kp
          n = n + 1
          sdata( 1, n ) = scale * sdata( 1, n )
          sdata( 2, n ) = scale * sdata( 2, n )
        end do
      end do
      return
      end
