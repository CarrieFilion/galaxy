      real*8 function dfmerr( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Osipkov-Merritt DF for a Jaffe sphere - as BT
c
c calling arguments
      real*8 E, Lz
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
      real*8 DawsnI, errfn
c
c local variables
      integer ifail, iuse
      real*8 arg, dfb, dfm, DFscle, Es, Escale, Ls, Lzscle
      real*8 q, t1, t2, t3, t4
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
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
c special formula for purely radial orbits
      if( dfm .eq. 0. )then
        dfmerr = 0
c purely radial orbits
        if( Lz .eq. 0. )then
          arg = sqrt( -Es )
          ifail = 0
          t4 = DawsnI( arg, ifail )
          arg = sqrt( -2. * Es )
          t2 = DawsnI( arg, ifail )
          dfmerr = ( sqrt( 2. ) * t4 - t2 ) / ( 2. * pi**3 )
        end if
      else
c normal form
        Ls = Lz * Lzscle
        q = -( Es + .5 * dfb * Ls * Ls )
        arg = sqrt( 2. * q )
        ifail = 0
        t1 = .5 * sqrt( pi ) * exp( arg * arg ) * errfn( arg, ifail )
        t2 = ( 1. + dfb ) * DawsnI( arg, ifail )
        arg = sqrt( q )
        t3 = .5 * sqrt( pi ) * exp( arg * arg ) * errfn( arg, ifail )
        t4 = ( 1. + .5 * dfb ) * DawsnI( arg, ifail )
        dfmerr = ( t1 + t2 - sqrt( 2. ) * ( t3 + t4 ) ) / ( 2 * pi**3 )
      end if
      dfmerr = DFscle * dfmerr
      return
      end
