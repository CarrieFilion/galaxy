      real*8 function dfnje( Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for sglint - the DF times the density of states function
c
c calling argument
      real*8 Lz
c
c common blocks
c
      include 'inc/params.f'
c
      common / cut / E
      real*8 E
c
      include 'inc/model.f'
c
c externals
      real*8 distfn, satab, tautab
c
c local variable
      include 'inc/pi.f'
c
c DFs for flattened spheroids
      if( sphrod( icmp ) )then
        dfnje = 4 * pi**2 * satab( E, Lz ) * distfn( E, Lz )
      else if( dist( icmp ) )then
c DFs for discs
        dfnje = 2. * pi * tautab( E, Lz ) * distfn( E, Lz )
c extra factor for spheres
        if( .not. disc( icmp ) )dfnje = 4. * pi * Lz * dfnje
      else
        call crash( 'DFNJE', 'Unknown DFTYPE' )
      end if
      return
      end
