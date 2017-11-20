      function glomin ( a, b, c, m, machep, e, t, f, x )

c*********************************************************************72
c
cc GLOMIN seeks a global minimum of a function F(X) in an interval [A,B].
c
c  Discussion:
c
c    This function assumes that F(X) is twice continuously differentiable
c    over [A,B] and that F''(X) <= M for all X in [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2008
c
c  Author:
c
c    Richard Brent
c    Modifications by John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the interval.
c    It must be the case that A < B.
c
c    Input, double precision C, an initial guess for the global
c    minimizer.  If no good guess is known, C = A or B is acceptable.
c
c    Input, double precision M, the bound on the second derivative.
c
c    Input, double precision MACHEP, an estimate for the relative machine
c    precision.
c
c    Input, double precision E, a positive tolerance, a bound for the
c    absolute error in the evaluation of F(X) for any X in [A,B].
c
c    Input, double precision T, a positive error tolerance.
c
c    Input, external double precision F, the name of a user-supplied
c    function, of the form "FUNCTION F ( X )", which evaluates the
c    function whose global minimum is being sought.
c
c    Output, double precision X, the estimated value of the abscissa
c    for which F attains its global minimum value in [A,B].
c
c    Output, double precision GLOMIN, the value F(X).
c
      implicit none

      double precision a
      double precision a0
      double precision a2
      double precision a3
      double precision b
      double precision c
      double precision d0
      double precision d1
      double precision d2
      double precision e
      double precision f
      double precision glomin
      double precision h
      integer k
      double precision m
      double precision m2
      double precision machep
      double precision p
      double precision q
      double precision qs
      double precision r
      double precision s
      double precision sc
      double precision t
      double precision x
      double precision y
      double precision y0
      double precision y1
      double precision y2
      double precision y3
      double precision yb
      double precision z0
      double precision z1
      double precision z2

      a0 = b
      x = a0
      a2 = a
      y0 = f ( b )
      yb = y0
      y2 = f ( a )
      y = y2

      if ( y0 .lt. y ) then
        y = y0
      else
        x = a
      end if

      if ( m .le. 0.0D+00 .or. a .ge. b ) go to 140

      m2 = 0.5D+00 * ( 1.0D+00 + 16.0D+00 * machep ) * m

      if ( sc .le. a .or. sc .ge. b ) then
        sc = 0.5D+00 * ( a + b )
      else
        sc = c
      end if

      y1 = f ( sc )
      k = 3
      d0 = a2 - sc
      h = 9.0D+00 / 11.0D+00

      if ( y1 .ge. y ) go to 30
      x = sc
      y = y1

30    continue

      d1 = a2 - a0
      d2 = sc - a0
      z2 = b - a2
      z0 = y2 - y1
      z1 = y2 - y0
      r = d1 * d1 * z0 - d0 * d0 * z1
      p = r
      qs = 2.0D+00 * ( d0 * z1 - d1 * z0 )
      q = qs
      if ( k .gt. 1000000 .and. y .lt. y2 ) go to 50

40    continue

      if ( q * ( r * ( yb - y2 ) + z2 * q * ( ( y2 - y ) + t ) ) .ge.
     &  z2 * m2 * r * ( z2 * q - r ) ) go to 50
      a3 = a2 + r / q
      y3 = f ( a3 )

      if ( y3 .lt. y ) then
        x = a3
        y = y3
      end if
c
c  Assume that 1611 * K does not overflow.
c
50    continue

      k = mod ( 1611 * k, 1048576 )
      q = 1.0D+00
      r = ( b - a ) * 0.00001D+00 * dble ( k )
      if ( r .lt. z2 ) go to 40
      r = m2 * d0 * d1 * d2
      s = sqrt ( ( ( y2 - y ) + t ) / m2 )
      h = 0.5D+00 * ( 1.0D+00 + h )
      p = h * ( p + 2.0D+00 * r * s )
      q = q + 0.5D+00 * qs
      r = - 0.5D+00 * ( d0 + ( z0 + 2.01D+00 * e ) / ( d0 * m2 ) )
      if ( r .ge. s .and. d0 .ge. 0.0D+00 ) go to 60
      r = a2 + s
      go to 70

60    continue

      r = a2 + r

70    continue

      if ( p * q .le. 0.0D+00 ) go to 80
      a3 = a2 + p / q
      go to 90

80    continue

      a3 = r

90    continue

      if ( a3 .lt. r ) a3 = r
      if ( a3 .lt. b ) go to 100
      a3 = b
      y3 = yb
      go to 110

100   continue

      y3 = f ( a3 )

110   continue

      if ( y3 .lt. y ) then
        x = a3
        y = y3
      end if

120   continue

      d0 = a3 - a2
      if ( a3 .le. r ) go to 130
      p = 2.0D+00 * ( y2 - y3 ) / ( m * d0 )
      if ( abs ( p ) .ge. ( 1.0D+00 + 9.0D+00 * machep ) * d0 ) 
     &  go to 130
      if ( 0.5D+00 * m2 * ( d0 * d0 + p * p ) .le.
     &  ( y2 - y ) + ( y3 - y ) + 2.0D+00 * t ) go to 130
      a3 = 0.5D+00 * ( a2 + a3 )
      h = 0.9D+00 * h
      go to 90

130   continue

      if ( a3 .ge. b ) go to 140
      a0 = sc
      sc = a2
      a2 = a3
      y0 = y1
      y1 = y2
      y2 = y3
      go to 30

140   continue

      glomin = y

      return
      end
