      real*8 function resrdf( r )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c returns kappa / m - omega - needed for resrad
c
c calling arguments
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      common / rdfres / am, rmn, ltheory
      logical ltheory
      real rmn
      real*8 am
c
c externals
      real*8 akappa, omegac
c
c local variables
      integer k
      real ak, dr, om, r1
c
      resrdf = 0
      if( r .ne. 0. )then
        if( ltheory )then
          resrdf = akappa( r ) / am - omegac( r )
        else
c get Omega and kappa from data array
          r1 = r * lscale
          k = nint( ( r1 - rmn ) / drfrqs )
          if( k .lt. mr - 2 )then
            dr = r1 - real( k ) * drfrqs
            om = ( 1. - dr ) * wres( k + 1 ) + dr * wres( k + 2 )
            k = k + mr
            ak = ( 1. - dr ) * wres( k + 1 ) + dr * wres( k + 2 )
          else
            dr = real( mr - 1 ) * drfrqs / lscale
            om = wres( mr - 1 ) * dr / r
            ak = wres( 2 * mr - 1 ) * dr / r
          end if
          resrdf = ak / am - om
        end if
      end if
      return
      end
