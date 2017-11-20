      function local_min ( a, b, eps, t, f, x )

c*********************************************************************72
c
cc LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
c
c  Discussion:
c
c    The method used is a combination of golden section search and
c    successive parabolic interpolation.  Convergence is never much slower
c    than that for a Fibonacci search.  If F has a continuous second
c    derivative which is positive at the minimum (which is not at A or
c    B), then convergence is superlinear, and usually of the order of
c    about 1.324....
c
c    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
c    F is never evaluated at two points closer than TOL.  
c
c    If F is a unimodal function and the computed values of F are always
c    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
c    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
c    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  
c
c    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
c    perhaps non-global, minimum to the same accuracy.
c
c    Thanks to Jonathan Eggleston for pointing out a correction to the 
c    golden section step, 01 July 2013.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2013
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
c
c    Input, double precision EPS, a positive relative error tolerance.
c    EPS should be no smaller than twice the relative machine precision,
c    and preferably not much less than the square root of the relative
c    machine precision.
c
c    Input, double precision T, a positive absolute error tolerance.
c
c    Input, external double precision F, the name of a user-supplied
c    function, of the form "FUNCTION F ( X )", which evaluates the
c    function whose local minimum is being sought.
c
c    Output, double precision X, the estimated value of an abscissa
c    for which F attains a local minimum value in [A,B].
c
c    Output, double precision LOCAL_MIN, the value F(X).
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      double precision e
      double precision eps
      double precision f
      double precision fu
      double precision fv
      double precision fw
      double precision fx
      double precision local_min
      double precision m
      double precision p
      double precision q
      double precision r
      double precision sa
      double precision sb
      double precision t
      double precision t2
      double precision tol
      double precision u
      double precision v
      double precision w
      double precision x
c
c  C is the square of the inverse of the golden ratio.
c
      c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

      sa = a
      sb = b
      x = sa + c * ( b - a )
      w = x
      v = w
      e = 0.0D+00
      fx = f ( x )
      fw = fx
      fv = fw

10    continue

      m = 0.5D+00 * ( sa + sb ) 
      tol = eps * abs ( x ) + t
      t2 = 2.0D+00 * tol
c
c  Check the stopping criterion.
c
      if ( abs ( x - m ) .le.  t2 - 0.5D+00 * ( sb - sa ) ) go to 190
      r = 0.0D+00
      q = r
      p = q
      if ( abs ( e ) .le. tol ) go to 40
c
c  Fit a parabola to the points.
c
      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )
      if ( q .le. 0.0D+00 ) go to 20
      p = - p
      go to 30

20    continue

      q = - q

30    continue

      r = e
      e = d

40    continue

      if ( abs ( p ) .ge. abs ( 0.5D+00 * q * r ) ) go to 60
      if ( p .le. q * ( sa - x ) .or. p .ge. q * ( sb - x ) ) go to 60
c
c  Take the parabolic interpolation step.
c
      d = p / q
      u = x + d
c
c  F must not be evaluated too close to A or B.
c
      if ( ( u - sa ) .ge. t2 .and. ( sb - u ) .ge. t2 ) go to 90
      if ( x .ge. m ) go to 50
      d = tol
      go to 90

50    continue

      d = - tol
      go to 90
c
c  A golden-section step.
c
60    continue

      if ( x .ge. m ) go to 70
      e = sb - x
      go to 80

70    continue

      e = sa - x

80    continue

      d = c * e
c
c  F must not be evaluated too close to X.
c
90    continue

      if ( abs ( d ) .lt. tol ) go to 100
      u = x + d
      go to 120

100   continue

      if ( d .le. 0.0D+00 ) go to 110
      u = x + tol
      go to 120

110   continue

      u = x - tol

120   continue

      fu = f ( u )
c
c  Update A, B, V, W, and X.
c
      if ( fu .gt. fx ) go to 150
      if ( u .ge. x ) go to 130
      sb = x
      go to 140

130   continue

      sa = x

140   continue

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 10

150   continue

      if ( u .ge. x ) go to 160
      sa = u
      go to 170

160   continue

      sb = u

170   continue

      if ( fu .gt. fw .and. w .ne. x ) go to 180
      v = w
      fv = fw
      w = u
      fw = fu
      go to 10

180   continue

      if ( fu .gt. fv .and. v .ne. x .and. v .ne. w ) go to 20
      v = u
      fv = fu
      go to 10

190   continue

      local_min = fx

      return
      end
