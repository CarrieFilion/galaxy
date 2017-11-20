      real*8 function dfDejo( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Dejonghe's distribution function for a Plummer sphere (MNRAS v224, p13, 1987)
c   input parameter q is simply related to beta = .5qr^2/(1+r^2) (eq 23)
c    q = 2 for maximal radial bias
c    0 < q < 2 radial bias
c    q = 0 for an isotropic model (a polytrope of index 5)
c    q < 0 for azimuthal bias
c    q = -infty for an Einstein sphere
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gam1o1, gam1o2, hyg2f1
c
c local variables
      integer ifail, iuse
      real*8 a, b, c, d, DFscle, Es, Escale, hyg, Ls, Lzscle, q, x
      include 'inc/pi.f'
      save iuse, DFscle, Escale, Lzscle
c
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
c allow for non-unit mass and/or scale for this component
        Escale = rscale( icmp ) / cmpmas( icmp )
        Lzscle = 1. / sqrt( rscale( icmp ) * cmpmas( icmp ) )
        DFscle = 1. / sqrt( rscale( icmp )**3 * cmpmas( icmp ) )
        iuse = 1
      end if
      Es = E * Escale
c
      q = dfcns( 1, icmp )
      if( q .eq. 0. )then
c isotropic model
        dfDejo = ( 3.d0 / ( 7.d0 * pi**3 ) ) * ( -2.d0 * Es )**3.5
      else if( q .lt. 1.9999d0 )then
c general anisotropic model
        Ls = Lz * Lzscle
        a = 0
        b = .5d0 * q
        c = 4.5d0 - q
        d = 1
        x = -Ls * Ls / ( 2.d0 * Es )
        ifail = 0
        if( x .le. 1.0d0 )then
          hyg = gam1o1( b, c ) * hyg2f1( b, 1.d0 - c, d, x )
        else
          x = 1.d0 / x
          hyg = gam1o2( b, d - b, b + c ) * x**b *
     +                                          hyg2f1( b, b, b + c, x )
        end if
        dfDejo = 1.5d0 * gam1o1( 6.d0 - q, .5d0 * q )
     +               * ( -Es )**( 3.5d0 - q ) * hyg / ( 2.d0 * pi )**2.5
      else
c maximum possible radial bias
        Ls = Lz * Lzscle
        x = -2.d0 * Es - Ls * Ls
        dfDejo = 0
        if( x .gt. 0.d0 )dfDejo = 6.d0 * x**1.5 / ( 2.d0 * pi )**3
      end if
      dfDejo = DFscle * dfDejo
      return
      end
