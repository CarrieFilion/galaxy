      real*8 function sigmav( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns dispersion of the azimuthal velocity components in a disk or sphere
c   at the requested radius
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
      real*8 sigmau, vmean, v2mean
c
c Kalnajs function
      if( ( cdft( icmp ) .eq. 'KALN' ) .and. kcold )then
        sigmav = 0
c self-similar MTZ disc - radius is a dummy argument
      else if( cdft( icmp ) .eq. 'ZANG' )then
        sigmav = 1. - vmean( r )**2
        sigmav = sqrt( sigmav )
c generic disc formula
      else if( disc( icmp ) )then
        sigmav = v2mean( r ) - vmean( r )**2
        sigmav = max( sigmav, 0.d0 )
        sigmav = sqrt( sigmav )
c Merritt function for Jaffe halo
      else if( cdft( icmp ) .eq. 'MERR' )then
        sigmav = sigmau( r )
        if( dfcns( 1, icmp ) .ge. 0. )sigmav =
     +                 sigmav / sqrt( 1. + ( r / dfcns( 1, icmp ) )**2 )
c Dejonghe function for Plummer halo - MNRAS v224, p13 (1987)
      else if( cdft( icmp ) .eq. 'DEJO' )then
        sigmav = sigmau( r ) *
     + sqrt( 1.d0 - .5d0 * dfcns( 1, icmp ) * r * r / ( 1.d0 + r * r ) )
c Polyachenko-Shukhman function for uniform sphere
      else if( cdft( icmp ) .eq. 'USPS' )then
        sigmav = 0.25 * cmpmas( icmp ) / rscale( icmp )
        sigmav = sqrt( sigmav )
      else
        call crash( 'SIGMAV', 'Unknown type of DF' )
      end if
      return
      end
