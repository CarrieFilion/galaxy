      real*8 function sigmau( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns dispersion of the radial velocity components in a disk or sphere
c   at the requested radius
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      external sigufn
      real*8 akappa, gsigma, isopsi, quad_Pat, qtoom, sgulia
c
c local variables
      integer i, ifail, n
      real*8 a, b, dfm, epsr, gsig, r1, r2, zero
      include 'inc/pi.f'
c
      dfm = dfcns( 1, icmp )
c vely dispersion is zero for a cold disc
      if( kcold )then
        sigmau = 0
c Kalnajs functions of the form: surface density = w**m * tau(y)
c  integral form - as eq.(23) of Kalnajs (1976) Ap. J. 205, 751
      else if( cdft( icmp ) .eq. 'KALN' )then
        epsr = 1.e-6
        ifail = 0
        zero = 0.
        r1 = .01
        r1 = max( r, r1 )
        gsig = gsigma( r1 )
        if( gsig .gt. 0.d0 )then
          sigmau = quad_Pat( sigufn, zero, r1, epsr, ifail )
          sigmau = -sigmau / ( gsig * r1**( dfm + 1. ) )
          sigmau = sqrt( sigmau )
        end if
c generalised Kalnajs function
      else if( cdft( icmp ) .eq. 'LIA ' )then
        sigmau = sgulia( r )
c Miyamoto's function for Kuz'min/Toomre disc
      else if( cdft( icmp ) .eq. 'MIYA' )then
        sigmau = ( 2. * dble( mmiy ) + 4. ) * sqrt( 1. + r * r )
        sigmau = 1. / sqrt( sigmau )
c Zang function in self-similar disc
      else if( cdft( icmp ) .eq. 'ZANG' )then
        sigmau = 1. / sqrt( 1. + dfm )
c Omega model
      else if( cdft( icmp ) .eq. 'OMEG' )then
        if( r .lt. 1.d0 )then
          sigmau = ( .75 * pi - omegak * omegak ) * ( 1. - r * r ) / 3.
          sigmau = sqrt( sigmau )
        else
          sigmau = 0
        end if
c Merritt function for Jaffe halo - formula (9) of MNRAS 214, 25p
      else if( cdft( icmp ) .eq. 'MERR' )then
        sigmau = .5
        if( r .gt. 0. )then
          b = dfcns( 2, icmp )
          r1 = 1. + r
          sigmau = (.5- 2.*r-(9.+1.5*b)*r*r-(6.+b)*r**3
     +          - r*r * r1*r1 * (6.+b) * log(r/r1)) / (1.+r*r*b)
          sigmau = sqrt( sigmau )
        end if
c Dejonghe function for Plummer halo - MNRAS 224, 13 (1987)
      else if( cdft( icmp ) .eq. 'DEJO' )then
        sigmau = ( 6.d+0 - dfm ) * sqrt( 1. + r * r )
        sigmau = 1. / sqrt( sigmau )
c disc function constructed by Lucy method
      else if( cdft( icmp ) .eq. 'LUCY' )then
        sigmau = qtoom( r ) * 3.3585 * gsigma( r ) / akappa( r )
c Neil's DF
      else if( cdft( icmp ) .eq. 'NEIL' )then
        sigmau = dfcns( 1, icmp ) * 3.358 * gsigma( r ) / akappa( r )
c Composite Omega model B0
      else if( cdft( icmp ) .eq. 'CPB0' )then
        if( r .lt. 1.d0 )then
          sigmau = ( 1. - r * r ) / pi
          sigmau = 4. * sqrt( sigmau ) / 3.
        else
          sigmau = 0
        end if
c spherical polytrope
      else if( cdft( icmp ) .eq. 'POLY' )then
        if( r .lt. 1.d0 )then
          sigmau = isopsi( r ) / ( dfcns( 3, icmp ) + 1. )
          sigmau = sqrt( sigmau )
        else
          sigmau = 0
        end if
c Evans's DFs for power law discs - preprint March 1993
      else if( cdft( icmp ) .eq. 'EVAN' )then
        if( r .gt. 0. )then
          sigmau = r**evbeta * ( 2. * evbeta + dfm + 1. )
          sigmau = 1. / sqrt( sigmau )
        else
          sigmau = 0
        end if
c Evans/Collett DFs for Rybicki discs - MN v264 p353
      else if( cdft( icmp ) .eq. 'EVCO' )then
        n = nint( dfm )
        r2 = r * r
        r1 = sqrt( 1. + r2 )
c Appendix C - the next few lines seem to be slightly incorrect
        if( n .ge. 0 )then
          sigmau = 1
          a = 2**n * r1 / ( 1. + r1 )**n
          do i = 0, n
            if( i .gt. 0 )a = a * r2 * dble( n + 1 - i )
     +                                 / ( ( 1. + r1 ) * dble( 2 * i ) )
            sigmau = sigmau + a / dble( n + i + 1 )
          end do
          sigmau = sigmau / ( 1. + r1 )
        else
c equation (2.8)
          sigmau = r1 * log( ( 1. + r1 ) / r1 )
        end if
        sigmau = sqrt( sigmau )
c singular isothermal sphere
      else if( cdft( icmp ) .eq. 'SISP' )then
        sigmau = sqrt( .5 * cmpmas( icmp ) )
c Polyachenko-Shukhman function for uniform sphere
      else if( cdft( icmp ) .eq. 'USPS' )then
        r1 = r / rscale( icmp )
        if( r1 .ge. 1.0 )then
          sigmau = 0.0
        else
          sigmau = .5 * cmpmas( icmp ) / rscale( icmp )
          sigmau = 0.5 * sigmau * ( 1. - r1 * r1 )
          sigmau = sqrt( sigmau )
        end if
      else
        call crash( 'SIGMAU', 'Unknown type of DF' )
      end if
      return
      end
