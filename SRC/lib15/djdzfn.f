      real*8 function djdzfn( t )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrand for dfdjdz - Dejonghe-deZeeuw DF for Kuz'min-Kutuzov spheroid
c   ref: Dejonghe & de Zeeuw - ApJ 333 p 104 (1988)
c
c calling argument
      real*8 t
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      common / djdz / eloc, eps, rz
      real*8 eloc, eps, rz
c
c local variables
      real*8 f1, f2, t1, xeps
c
      t1 = sqrt( 1. - t * t )
      xeps = 2. * eloc * t * t1 / ( 1. + eps * t * rz )
      f1 = t1 * t1 / ( 1. - 2. * eloc * t * t1 + eps * t * rz )**5
      f2 = ( 3. + 4. * xeps - xeps * xeps ) * ( 1. - xeps ) * t1 * t1 +
     +     12. * t * t
      djdzfn = f1 * f2
      return
      end
