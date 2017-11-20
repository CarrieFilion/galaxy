      real*8 function rannor( mean, sigma )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns a random number from a normal distribution centred on mean with
c   dispersion sigma
c has same calling arguments and functionality as NAg routine g05ddf
c
c calling arguments
      real*8 mean, sigma
c
c external
      real*8 ranuni
c
c local variables
      integer j
      real*8 a( 2 ), x, y
      save a, j
      include 'inc/pi.f'
      data j / 2 /
c
      j = j + 1
      if( j .gt. 2 )then
c von Neumann's algorithm
        x = 2. * pi * ranuni( 0. )
        y = ranuni( 0. )
        y = -2. * log( y )
        y = sqrt( y )
        a( 1 ) = sin( x ) * y
        a( 2 ) = cos( x ) * y
        j = 1
      end if
      rannor = mean + sigma * a( j )
      return
      end
