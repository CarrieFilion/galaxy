      real*8 function isrhop( psi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c evaluates density as function of the relative potential psi for spherical
c   isotropic models
c
c calling argument
      real*8 psi
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c external - error function
      real*8 errfn
c
c local variables
      integer ifail
      real*8 arg, t
      include 'inc/pi.f'
c
      isrhop = 0
c inside the tidal radius
      if( psi .gt. 0.d0 )then
c King model - expression (4-131) from Binney & Tremaine
        if( ctype( icmp ) .eq. 'KING' )then
          arg = psi / sigmak**2
          ifail = 0
          t = sqrt( arg )
          isrhop = exp( arg ) * errfn( t, ifail ) -
     +                    2.d0 * t * ( 1.d0 + arg / 1.5d0 ) / sqrt( pi )
          isrhop = rhofac * isrhop
c polytrope - expression (4-107a) from Binney & Tremaine
        else if( ctype( icmp ) .eq. 'POLY' )then
          isrhop = rhofac * psi**dfcns( 3, icmp )
        else
          call crash( 'ISRHOP', 'Unrecognized halo' )
        end if
      end if
      return
      end
