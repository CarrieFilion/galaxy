      real*8 function phisph( R, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns potential at point (R,z) for an axisymmetric mass distribution
c
c calling arguments
      real*8 R, z
c
c common blocksc
c
      include 'inc/params.f'

      include 'inc/model.f'
c
c externals
      real*8 phikks, phispl
c
c local variables
      real*8 t
c
c spline fit potential from DFITER
      if( ctype( icmp ) .eq. 'DFIT' )then
        phisph = phispl( R, z )
c Kuz'min-Kutuzov oblate sphere
      else if( ctype( icmp ) .eq. 'KKSP' )then
        phisph = phikks( R, z )
c Miyamoto-Nagai model - eq. (2-50a) of BT
      else if( ctype( icmp ) .eq. 'MYNG' )then
        t = sqrt( z**2 + dfcns( 3, icmp )**2 )
        phisph = -1. / sqrt( R**2 + ( 1. + t )**2 )
c none of the above
      else
        call crash( 'PHISPH', 'Unknown type of spheroid' )
      end if
      return
      end
