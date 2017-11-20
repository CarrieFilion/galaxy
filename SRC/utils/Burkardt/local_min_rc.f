      subroutine local_min_rc ( a, b, arg, status, value )

c*********************************************************************72
c
cc LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
c
c  Discussion:
c
c    This routine seeks an approximation to the point where a function
c    F attains a minimum on the interval (A,B).
c
c    The method used is a combination of golden section search and
c    successive parabolic interpolation.  Convergence is never much
c    slower than that for a Fibonacci search.  If F has a continuous
c    second derivative which is positive at the minimum (which is not
c    at A or B), then convergence is superlinear, and usually of the
c    order of about 1.324...
c
c    The routine is a revised version of the Brent local minimization 
c    algorithm, using reverse communication.
c
c    It is worth stating explicitly that this routine will NOT be
c    able to detect a minimizer that occurs at either initial endpoint
c    A or B.  If this is a concern to the user, then the user must
c    either ensure that the initial interval is larger, or to check
c    the function value at the returned minimizer against the values
c    at either endpoint.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    Algorithms for Minimization Without Derivatives,
c    Dover, 2002,
c    ISBN: 0-486-41998-3,
c    LC: QA402.5.B74.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters
c
c    Input/output, double precision A, B.  On input, the left and right
c    endpoints of the initial interval.  On output, the lower and upper
c    bounds for an interval containing the minimizer.  It is required
c    that A < B.
c
c    Output, double precision ARG, the currently considered point.  The user
c    does not need to initialize this value.  On return with STATUS positive,
c    the user is requested to evaluate the function at ARG, and return
c    the value in VALUE.  On return with STATUS zero, ARG is the routine's
c    estimate for the function minimizer.
c
c    Input/output, integer STATUS, used to communicate between the user
c    and the routine.  The user only sets STATUS to zero on the first call,
c    to indicate that this is a startup call.  The routine returns STATUS
c    positive to request that the function be evaluated at ARG, or returns
c    STATUS as 0, to indicate that the iteration is complete and that
c    ARG is the estimated minimizer.
c
c    Input, double precision VALUE, the function value at ARG, as requested
c    by the routine on the previous call.
c
c  Local parameters:
c
c    C is the squared inverse of the golden ratio.
c
c    EPS is the square root of the relative machine precision.
c
      implicit none

      double precision a
      double precision arg
      double precision b
      double precision c
      double precision d
      double precision e
      double precision eps
      double precision fu
      double precision fv
      double precision fw
      double precision fx
      double precision midpoint
      double precision p
      double precision q
      double precision r
      double precision r8_epsilon
      integer status
      double precision tol
      double precision tol1
      double precision tol2
      double precision u
      double precision v
      double precision value
      double precision w
      double precision x

      save c
      save d
      save e
      save eps
      save fu
      save fv
      save fw
      save fx
      save midpoint
      save p
      save q
      save r
      save tol
      save tol1
      save tol2
      save u
      save v
      save w
      save x
c
c  STATUS (INPUT) = 0, startup.
c
      if ( status .eq. 0 ) then

        if ( b .le. a ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LOCAL_MIN_RC - Fatal error!'
          write ( *, '(a)' ) '  A < B is required, but'
          write ( *, '(a,g14.6)' ) '  A = ', a
          write ( *, '(a,g14.6)' ) '  B = ', b
          status = -1
          stop
        end if

        c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

        eps = sqrt ( r8_epsilon ( ) )
        tol = r8_epsilon ( )

        v = a + c * ( b - a )
        w = v
        x = v
        e = 0.0D+00

        status = 1
        arg = x

        return
c
c  STATUS (INPUT) = 1, return with initial function value of FX.
c
      else if ( status .eq. 1 ) then

        fx = value
        fv = fx
        fw = fx
c
c  STATUS (INPUT) = 2 or more, update the data.
c
      else if ( 2 .le. status ) then

        fu = value

        if ( fu .le. fx ) then

          if ( x .le. u ) then
            a = x
          else
            b = x
          end if

          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu

        else

          if ( u .lt. x ) then
            a = u
          else
            b = u
          end if

          if ( fu .le. fw .or. w .eq. x ) then
            v = w
            fv = fw
            w = u
            fw = fu
          else if ( fu .le. fv .or. v .eq. x .or. v .eq. w ) then
            v = u
            fv = fu
          end if

        end if

      end if
c
c  Take the next step.
c
      midpoint = 0.5D+00 * ( a + b )
      tol1 = eps * abs ( x ) + tol / 3.0D+00
      tol2 = 2.0D+00 * tol1
c
c  If the stopping criterion is satisfied, we can exit.
c
      if ( abs ( x - midpoint ) .le. 
     &   ( tol2 - 0.5D+00 * ( b - a ) ) ) then
        status = 0
        return
      end if
c
c  Is golden-section necessary?
c
      if ( abs ( e ) .le. tol1 ) then
        if ( midpoint .le. x ) then
          e = a - x
        else
          e = b - x
        end if

        d = c * e
c
c  Consider fitting a parabola.
c
      else

        r = ( x - w ) * ( fx - fv )
        q = ( x - v ) * ( fx - fw )
        p = ( x - v ) * q - ( x - w ) * r
        q = 2.0D+00 * ( q - r )
        if ( 0.0D+00 .le. q ) then
          p = - p
        end if
        q = abs ( q )
        r = e
        e = d
c
c  Choose a golden-section step if the parabola is not advised.
c
        if ( 
     &    ( abs ( 0.5D+00 * q * r ) .le. abs ( p ) ) .or. 
     &    ( p .le. q * ( a - x ) ) .or. 
     &    ( q * ( b - x ) .le. p ) ) then

          if ( midpoint .le. x ) then
            e = a - x
          else
            e = b - x
          end if

          d = c * e
c
c  Choose a parabolic interpolation step.
c
        else

          d = p / q
          u = x + d

          if ( ( u - a ) .lt. tol2 ) then
            d = sign ( tol1, midpoint - x )
          end if

          if ( ( b - u ) .lt. tol2 ) then
            d = sign ( tol1, midpoint - x )
          end if

        end if

      end if
c
c  F must not be evaluated too close to X.
c
      if ( tol1 .le. abs ( d ) ) then
        u = x + d
      end if

      if ( abs ( d ) .lt. tol1 ) then
        u = x + sign ( tol1, d )
      end if
c
c  Request value of F(U).
c
      arg = u
      status = status + 1

      return
      end
