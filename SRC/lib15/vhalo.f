      real*8 function vhalo( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns circular velocity due to the halo only
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
c externals
      real*8 frtot, gmassh, vbss, vcirc, vdisc
c
c local variables
      integer ip, jcmp
      real*8 dmass, haloc, hmass, hrad, r1, vt, vtot2
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'VHALO',
     +                                        'Nonsense value of icmp' )
      if( disc( icmp ) )then
        call crash( 'VHALO',
     +                     'Halo function called for a disc component' )
      else
        hmass = cmpmas( icmp )
        haloc = dfcns( 3, icmp )
        hrad = rscale( icmp )
      end if
      vhalo = 0
c no halo
      if( ( r .eq. 0. ) .or. ( disc( icmp ) ) )then
        vhalo = 0.
c frozen disc model
      else if( ctype( icmp ) .eq. 'ASDI' )then
        dmass = 0
        jcmp = icmp
        do ip = 1, ncmp
          dmass = dmass + cmpmas( ip )
          if( disc( ip ) )icmp = ip
        end do
        vhalo = sqrt( hmass / dmass ) * vdisc( r )
        icmp = jcmp
c uniform sphere
      else if( ctype( icmp ) .eq. 'UNIS' )then
        r1 = abs( r ) / hrad
        if( r1 .le. 1. )then
          vhalo = r1 * sqrt( hmass / hrad )
        else
          vhalo = sqrt( hmass / abs( r ) )
        end if
c Plummer sphere
      else if( ctype( icmp ) .eq. 'PLUM' )then
        vhalo = r * sqrt( hmass ) / ( hrad * hrad + r * r )**.75
c Fall/Efstathiou model
      else if( ctype( icmp ) .eq. 'FE  ' )then
        r1 = r / hrad
        vtot2 = haloc * haloc * r1 * r1 / ( 1. + r1 * r1 )
        jcmp = icmp
        do ip = 1, ncmp
          if( disc( ip ) )icmp = ip
        end do
        vt = vtot2 - vdisc( r )**2
        icmp = jcmp
        if( vt .gt. 0. )then
          vhalo = sqrt( vt )
        else
          vhalo = -sqrt( -vt )
        end if
c cored isothermal halo
      else if( ctype( icmp ) .eq. 'ISOT' )then
        r1 = abs( r ) / hrad
c        vhalo = haloc * sqrt( 1. - atan( r1 ) / r1 )
        vhalo = cmpmas( icmp ) * r1**2 / ( hrad * ( 1. + r1**2 ) )
        vhalo = sqrt( vhalo )
c Bahcall/Schmidt/Soniera model
      else if( ctype( icmp ) .eq. 'BSS ' )then
        vt = vbss( abs( r ) )
        vhalo = sqrt( vt**2 - vdisc( r )**2 )
c NGC 3198 model
      else if( ctype( icmp ) .eq. '3198' )then
        vhalo = vcirc( r )**2 - vdisc( r )**2
        vhalo = sqrt( vhalo )
c Jaffe model
      else if( ctype( icmp ) .eq. 'JAFF' )then
        vhalo = sqrt( hmass / ( hrad + r ) )
c Kepler
      else if( ctype( icmp ) .eq. 'KEPL' )then
        r1 = max( abs( r ), 1.d - 12 )
        vhalo = sqrt( hmass / r1 )
c DFIT potential
      else if( ctype( icmp ) .eq. 'DFIT' )then
        vtot2 = -r * frtot( r )
        jcmp = icmp
        do icmp = 1, ncmp
          if( disc( icmp ) )vtot2 = vtot2 - vdisc( r )**2
        end do
        icmp = jcmp
        vhalo = max( vhalo, vtot2 )
        vhalo = sqrt( vhalo )
c the TEDS model
      else if( ctype( icmp ) .eq. 'TEDS' )then
        r1 = haloc * hrad
        vhalo = .5 / ( hrad * hrad + r * r )**1.5 +
     +          .09 / ( r1 * r1 + r * r )**1.5
        vhalo = r * sqrt( hmass * vhalo / .59 )
c Kuz'min-Kutuzov oblate spheroid - value in mid-plane
      else if( ctype( icmp ) .eq. 'KKSP' )then
        r1 = sqrt( hrad * hrad + r * r )
        vhalo = sqrt( hmass ) * r /
     +                            ( sqrt( r1 ) * ( hrad * haloc + r1 ) )
c Hernquist model
      else if( ctype( icmp ) .eq. 'HERN' )then
        r1 = abs( r )
        vhalo = sqrt( hmass * r1 ) / ( hrad + r1 )
c use Plummer sphere expression as default for an unknown model
      else if( ctype( icmp ) .eq. 'UNKN' )then
        vhalo = r * sqrt( hmass ) / ( hrad * hrad + r * r )**.75
c generalized Aguilar-Merritt model
      else if( ctype( icmp ) .eq. 'AGME' )then
        r1 = r / hrad
        if( r1 .lt. 1. )then
          vhalo = sqrt( hmass * r1**( 2. + haloc ) )
        else
          vhalo = sqrt( hmass / r1 )
        end if
c Henon's spherical isochrone
      else if( ctype( icmp ) .eq. 'ISOC' )then
        r1 = sqrt( hrad * hrad + r * r )
        vhalo = sqrt( hmass ) * r / ( sqrt( r1 ) * ( hrad + r1 ) )
c      else if( ctype( icmp ) .eq. 'NFW ' )then
c NFW profile
c        r1 = abs( r ) / hrad + 1.d-10
c        vhalo = log( 1. + r1 ) / r1**2 - 1. / ( r1 * ( 1. + r1 ) )
c        vhalo = sqrt( hmass * abs( r ) / hrad**2 * abs( vhalo ) )
c        vhalo = sign( vhalo, r )
c singular isothermal sphere
      else if( ctype( icmp ) .eq. 'SISP' )then
        vhalo = sqrt( hmass )
c all other models
      else
        vhalo = sqrt( gmassh( abs( r ) ) / abs( r ) )
      end if
      vhalo = sign( vhalo, r )
      return
      end
