      real*8 function gofxkt( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Kalnajs eq. (34) for Kuz'min/Toomre disc - power series in x**2
c
c calling argument
      real*8 x
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      integer j
      real*8 a, dfm, pf, term, xf
      include 'inc/pi.f'
c
      dfm = dfcns( 1, icmp )
c j = 0 term is unity
      j = 0
      gofxkt = 1.
      xf = 1.
      pf = 1
      term = 1
c convergence test - function and terms are all positive
      do while ( term .gt. 1.d-8 * gofxkt )
        a = j
        j = j + 1
        pf = pf * ( ( .5 * dfm + .5 + a ) / ( .5 + a ) )
        pf = pf * ( ( .5 * dfm + 1. + a ) / ( dfm + a ) )
        pf = pf * ( ( .5 * dfm - 1.5 + a ) / ( 1. + a ) )
        xf = xf * x * x
        term = pf * xf
        gofxkt = gofxkt + term
      end do
c constant factor
      gofxkt = gofxkt * dfm / ( 2. * pi**2 )
      return
      end
