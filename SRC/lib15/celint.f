      real*8 function celint( kc, pp, aa, bb )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
C Complete Elliptic Integral - algorithm and notation as given in Numerical
C         Recipes pp183-188
c
c calling arguments
      real*8 aa, bb, kc, pp
c
c local variables
      real*8 a, b, ca, e, em, f, g, p, q, qc
c accuracy = ca**2
      parameter ( ca = .000001 )
      include 'inc/pi.f'
c
      qc = abs( kc )
      a = aa
      b = bb
      p = pp
      e = qc
      em = 1
      if( p .gt. 0. )then
        p = sqrt( p )
        b = b / p
      else
        f = qc * qc
        q = 1. - f
        g = 1. - p
        f = f - p
        q = q * ( b - a * p )
        p = sqrt( f / g )
        a = ( a - b ) / g
        b = -q / ( g * g * p ) + a * p
      end if
    1 f = a
      a = a + b / p
      g = e / p
      b = b + f * g
      b = b + b
      p = g + p
      g = em
      em = qc + em
      if( abs( g - qc ) .gt. g * ca )then
        qc = sqrt( e )
        qc = qc + qc
        e = qc * em
        go to 1
      end if
      celint = .5 * pi * ( b + a * em ) / ( em * ( em + p ) )
      return
      end
