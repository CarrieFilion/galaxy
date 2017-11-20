      real*8 function gmassd( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns disc mass interior to r
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gmfund, gmtabd
c
c local arrays
      real*8, allocatable :: absc( : )
      real*8, allocatable :: wt( : )
c
c local variables
      integer i, ipar, npts
      real*8 cc, dmass, x, x1, x2
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'GMASSD',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        dmass = cmpmas( icmp )
        cc = dfcns( 3, icmp )
      else
        call crash( 'GMASSD', 'Called for a non-disc component' )
      end if
      x = r / rscale( icmp )
c
      if( r .le. rhole )then
        gmassd = 0
c Kuz'min/Toomre disc
      else if( ( ctype( icmp ) .eq. 'KT  ' ) .or.
     +         ( ctype( icmp ) .eq. 'SOFT' ) )then
        gmassd = dmass * ( 1. / sqrt( 1. + rhole * rhole )
     +                    - 1. / sqrt( 1. + x * x ) )
c Maclaurin/Freeman/Kalnajs disc
      else if( ctype( icmp ) .eq. 'MFK ' )then
        gmassd = ( 1. - rhole * rhole )**1.5
        if( x .lt. 1. )then
          gmassd = gmassd - ( 1. - x * x )**1.5
        end if
        gmassd = dmass * gmassd
c Mestel/Toomre/Zang disc
      else if( ctype( icmp ) .eq. 'MTZ ' )then
        gmassd = dmass * ( x - rhole )
c exponential disc or Sanders model
      else if( ( ctype( icmp ) .eq. 'EXP ' ) .or.
     +         ( ctype( icmp ) .eq. 'SAND' ) )then
        if( rhole .gt. 0. )then
          gmassd = ( 1. + rhole ) * exp( -rhole )
        else
          gmassd = 1
        end if
        gmassd = gmassd - ( 1. + x ) * exp( -x )
        gmassd = dmass * gmassd
c Hunter discs
      else if( ctype( icmp ) .eq. 'HUNT' )then
        ipar = impar( icmp )
        if( rhole .gt. 0. )then
          gmassd = ( 1. - rhole * rhole )**( real( ipar ) + .5 )
        else
          gmassd = 1
        end if
        if( x .lt. 1. )then
          gmassd = gmassd - ( 1. - x * x )**( real( ipar ) + .5 )
        end if
        gmassd = dmass * gmassd
c power law disc
      else if( ctype( icmp ) .eq. 'POWE' )then
        x1 = 2.d0 * cc + 1.d0
        gmassd = 2.d0 * pi * evsigma * ( x**x1 - rhole**x1 ) / x1
        gmassd = dmass * gmassd
c "cosine" disc - cc = x / a
      else if( ctype( icmp ) .eq. 'COSI' )then
        x1 = ( x - 1. ) * cc
        if( x1 .lt. -1. )then
          gmassd = 0
        else if( x1 .ge. 1. )then
          gmassd = dmass
        else
          x2 = 1. - x1 * x1
          gmassd = .5 * x1 * sqrt( x2 )
     +           + .5 * asin( x1 ) + .25 * pi - x2**1.5 / ( 3. * cc )
          gmassd = 2. * dmass * gmassd / pi
        end if
c Gaussian disc
      else if( ctype( icmp ) .eq. 'GAUS' )then
        if( rhole .gt. 0. )then
          gmassd = exp( -.5 * rhole**2 )
        else
          gmassd = 1
        end if
        gmassd = gmassd - exp( -.5 * x**2 )
        gmassd = dmass * gmassd
c finite Mestel disc
      else if( ctype( icmp ) .eq. 'FMES' )then
        if( x .lt. 1. )then
          gmassd = x * ( .5 * pi - asin( x ) ) + 1. - sqrt( 1. - x**2 )
        else
          gmassd = 1
        end if
        gmassd = dmass * gmassd
c Rybicki disc
      else if( ctype( icmp ) .eq. 'RYBI' )then
        gmassd = dmass * ( sqrt( 1. + x * x ) - 1. )
c composite model
      else if( ctype( icmp ) .eq. 'COMP' )then
c exponential part
        gmassd = 1. - ( 1. + x ) * exp( -x )
c Gaussian part
        x1 = 3 * x
        gmassd = gmassd + .1 * ( 1. - exp( -x1**2 ) )
        gmassd = dmass * gmassd
c tabulated disk
      else if( ctype( icmp ) .eq. 'MTAB' )then
        gmassd = gmtabd( r )
c Donner-Thomasson disk A&A v290 p785
      else if( ctype( icmp ) .eq. 'DOTH' )then
        gmassd = 1 + ( 1. + 2. * x ) * exp( -2. * x ) / 3. -
     +           4. * ( 1. + x ) * exp( -x ) / 3.
        gmassd = dmass * gmassd
      else
c arbitrary disc: non-adaptive Gauss-Legendre quadrature
        x1 = rhole
        x2 = min( r, rmax )
        npts = 32
        allocate ( absc( npts ) )
        allocate ( wt( npts ) )
        call GLqtab( x1, x2, npts, wt, absc )
        gmassd = 0
        do i = 1, npts
          x1 = absc( i )
          gmassd = gmassd + wt( i ) * gmfund( x1 )
        end do
      end if
      return
      end
