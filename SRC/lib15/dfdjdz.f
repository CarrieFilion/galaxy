      real*8 function dfdjdz( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns value of the Dejonghe-deZeeuw DF for a given E & Lz in the
c   potential of the Kuz'min-Kutuzov spheroidal isochrone
c   ref: eq (4.28) of Dejonghe & de Zeeuw - ApJ 333 p 104 (1988)
c   They use units s.t. phi(0) = 1, so I have to divide E by phi(0) and
c   Lz by sqrt( a + c ) [one factor for the distance scale less a square
c   root factor for the velocity scale, since phi(0) = -1 / ( a + c ) ]
c
c calling arguments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      common / djdz / Eloc, eps, rz
      real*8 Eloc, eps, rz
c
c externals
      external djdzfn
      real*8 quad_Pat
c
c local variables
      integer ifail
      real*8 a, b, epsr, fint1, fint2, haloc
      include 'inc/pi.f'
c
      haloc = dfcns( 3, icmp )
c store values for integrand - Eloc = their a * E
      Eloc = E / ( Emine * ( 1. + haloc ) )
c normalizing factors cancel in this expression
      rz = sqrt( -2. * ( 1. - haloc**2 ) * E * Lz * Lz )
c compute two integrals
      a = 0
      b = 1
      epsr = 1.e-6
      ifail = 0
      eps = -1.
      fint1 = quad_Pat( djdzfn, a, b, epsr, ifail )
      eps = 1.
      fint2 = quad_Pat( djdzfn, a, b, epsr, ifail )
c final expression
      dfdjdz = haloc**2 * ( E / Emine )**2.5 * ( fint1 + fint2 ) /
     +                           ( ( 1. + haloc ) * sqrt( 2. ) * pi )**3
c M must be divided by ( a + c ) and the phase space element is
c   multiplied by ( a + c )**3/2
      dfdjdz = dfdjdz * sqrt( 1. + haloc )
      return
      end
