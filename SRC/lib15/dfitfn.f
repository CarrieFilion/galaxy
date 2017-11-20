      real*8 function dfitfn( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c type of DF for halo plus fixed external potential - isotropic for now
c   may use NAG routines ERRFN and DAWSNI
c
c calling arguments
      real*8 E
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
      integer ifail
      real*8 arg, dfb, dfm, haloc, q, q2, tiny, t1, t2, t3, t4
      parameter ( tiny = 1.d-10 )
      include 'inc/pi.f'
c
      if( cdft( icmp ) .ne. 'DFIT' )call crash( 'DFITFN',
     +                                                 'Wrong DF type' )
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
      haloc = dfcns( 3, icmp )
c Polytropic form for DF iteration
      if( impar( icmp ) .eq. 1 )then
        if( dfm .eq. 0.d0 )then
          dfitfn = dfb
        else
          dfitfn = max( haloc - E, 0.d0 )
          dfitfn = dfb * dfitfn**dfm
        end if
      else if( impar( icmp ) .eq. 2 )then
c lowered isothermal (King model) form for DF iteration
        dfitfn = max( haloc - E, 0.d0 )
        dfitfn = dfb * ( exp( dfitfn / sigmai2 ) - 1.d0 )
      else if( impar( icmp ) .eq. 3 )then
c isotropic Jaffe model form for DF iteration
        dfitfn = 0
        if( E .lt. Phimax )then
          arg = sqrt( 2. * ( Phimax - E ) )
          ifail = 0
          t1 = 
     +         .5 * sqrt( pi ) * exp( arg * arg ) * errfn( arg, ifail )
          t2 = DawsnI( arg, ifail )
          arg = sqrt( Phimax - E )
          t3 =
     +         .5 * sqrt( pi ) * exp( arg * arg ) * errfn( arg, ifail )
          t4 = DawsnI( arg, ifail )
          dfitfn =
     +            ( t1 + t2 - sqrt( 2. ) * ( t3 + t4 ) ) / ( 2 * pi**3 )
          dfitfn = dfb * dfitfn
        end if
      else if( impar( icmp ) .eq. 4 )then
c isotropic Hernquist model for DF iteration
        q2 = E / Emine
        if( q2 .gt. 0. )then
          q = sqrt( q2 )
          dfitfn = ( 3. * asin( q ) + q * sqrt( 1. - q2 ) *
     +             ( 1. - 2. * q2 ) * ( 8. * q2 * ( q2 - 1. ) - 3. ) ) /
     +             ( 8. * sqrt( 2. ) * pi**3 * ( 1. - q2 + tiny )**2.5 )
          dfitfn = dfb * dfitfn
        else
          dfitfn = 0
        end if
      else
c none of the above
        print *, 'impar( icmp ) =', impar( icmp )
        call crash( 'DFITFN', 'Unknown function type' )
      end if
      return
      end
