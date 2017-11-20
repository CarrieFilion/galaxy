      real*8 function vmagris( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns value of mean orbital velocity for Kalnajs DF in a possibly
c   softened KT disk, as given by eq. (33)
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      integer k, l, m
      logical ksum, lsum
      real*8 c2, fk, fk2, fl, fl2, r2, sumk, suml, t, tk, tl, tl1, w, yy
      real*8 dfm
      include 'inc/pi.f'
c
      dfm = dfcns( 1, icmp )
      m = dfm + .5
      c2 = rstar * rstar
      r2 = r * r
      yy = r2 / ( c2 + r2 )
      w = rstar / sqrt( c2 + r2 )
      vmagris = 0.
      k = 0
      sumk = 0.
      fk = k
      fk2 = 2 * k
      tk = 1.
      tl1 = 2. / pi
      do l = 1, m
        fl = l
        tl1 = tl1 * fl / ( fl + .5 )
      end do
c sum over k
      ksum = .true.
      do while ( ksum )
        l = 0
        tl = tl1
        suml = tl
c sum over l
        lsum = .true.
        do while ( lsum )
          l = l + 1
          fl = l
          fl2 = 2 * l
          tl = tl * ( ( .5 * dfm - 2.5 + fl ) / fl )
          tl = tl * ( ( dfm + fk2 + fl2 ) / ( fk + fl - .5 ) )
          tl = tl * ( ( fk + fl ) / ( dfm + fk2 + fl2 + .5 ) )
          tl = tl * ( ( dfm + fk2 + fl2 - 1. ) /
     +                              ( dfm + fk2 + fl2 - .5 ) )
          tl = tl * yy
          suml = suml + tl
          lsum = tl .gt. 0.1d-3 * suml
        end do
        t = tk * suml
        sumk = sumk + t
        vmagris = sumk - .5 * t
        ksum = abs( t ) .gt. abs( 0.2d-3 * vmagris )
c next k
        if( ksum )then
          k = k + 1
          fk = k
          fk2 = 2 * k
c tk will be zero here for the hard potential
          tk = -tk * ( c2 - 1. ) * yy * ( fk + .5 ) / fk
c reset first term of l sum
          tl1 = tl1 * ( ( dfm + fk2 ) / ( fk - .5 ) )
          tl1 = tl1 * ( ( dfm + fk2 - 1. ) / ( dfm + fk2 - .5 ) )
          tl1 = tl1 * ( fk / ( dfm + fk2 + .5 ) )
        end if
      end do
      vmagris = vmagris * w**dfm
      vmagris = vmagris * ( 1. + r2 )**1.5 * sqrt( 2. * w / rstar )
      return
      end
