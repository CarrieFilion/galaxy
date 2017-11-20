      real*8 function confhg( a, b, x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Kummer's confluent hypergeometric function M( a, b, x )
c   or equivalently f( a; b; x )
c
c calling arguments
      real*8 a, b, x
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      common / confhn / aa, bb, z
      real*8 aa, bb, z
c
c externals
      external conffn
      real*8 Gammaf, quad_Pat
c
c local variables
      integer ifail
      real*8 epsr, one, rel, zero
      parameter ( zero = 0, one = 1 )
c
      aa = a
      bb = b
      z = x
c integral formula 13.2.1, p505 of Abramowitz & Stegun
      epsr = 1.d-8
      ifail = 0
      confhg = quad_Pat( conffn, zero, one, epsr, ifail )
c Gamma function
      confhg = confhg * Gammaf( b, ifail ) /
     +                   ( Gammaf( a, ifail ) * Gammaf( b - a, ifail ) )
      return
      end

      real*8 function conffn( t )
      real*8 t
c
c common block
c
      common / confhn / aa, bb, z
      real*8 aa, bb, z
c
      conffn = exp( z * t ) * t**( aa - 1 ) *
     +         ( 1. - t )**( bb - aa - 1 )
      return
      end
