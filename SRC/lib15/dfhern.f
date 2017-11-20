      real*8 function dfhern( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Hernquist's isotropic DF for his bulge/halo model (ApJ 1990 v356 p359)
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
c local variables
      integer iuse
      real*8 dfm, DFscle, Es, Escale, Ls, Lzscle, q, q2, tiny
      parameter ( tiny = 1.d-10 )
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
      Ls = Lz * Lzscle
c
      dfm = dfcns( 1, icmp )
c isotropic case
      if( dfm .le. 0. )then
        q2 = -Es
      else
c the anisotropy radius is stored as dfm
        q2 = -Es - Ls**2 / ( 2. * dfm**2 )
        call crash( 'DFHERN',
     +    'Anisotropic Hernquist models require a different potential' )
      end if
      dfhern = 0
      if( q2 .gt. 0. )then
        q = sqrt( q2 )
        dfhern = ( 3. * asin( q ) + q * sqrt( 1. - q2 ) *
     +             ( 1. - 2. * q2 ) * ( 8. * q2 * ( q2 - 1. ) - 3. ) ) /
     +             ( 8. * sqrt( 2. ) * pi**3 * ( 1. - q2 + tiny )**2.5 )
c        if( dfm .gt. 0. )dfhern = dfhern +
c     +            q * ( 1. - 2. * q2 ) / ( sqrt( 2. ) * pi**3 * dfm**2 )
      end if
      dfhern = DFscle * dfhern
      return
      end
