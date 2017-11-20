      subroutine zero_rc ( a, b, t, arg, status, value )

c*********************************************************************72
c
cc ZERO_RC seeks the root of a function F(X) using reverse communication.
c
c  Discussion:
c
c    The interval [A,B] must be a change of sign interval for F.
c    That is, F(A) and F(B) must be of opposite signs.  Then
c    assuming that F is continuous implies the existence of at least
c    one value C between A and B for which F(C) = 0.
c
c    The location of the zero is determined to within an accuracy
c    of 6 * MACHEPS * abs ( C ) + 2 * T.
c
c    The routine is a revised version of the Brent zero finder 
c    algorithm, using reverse communication.
c
c    Thanks to Thomas Secretin for pointing out a transcription error in the
c    setting of the value of P, 11 February 2013.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 February 2013
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
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the change of sign interval.
c
c    Input, double precision T, a positive error tolerance.
c
c    Output, double precision ARG, the currently considered point.  The user
c    does not need to initialize this value.  On return with STATUS positive,
c    the user is requested to evaluate the function at ARG, and return
c    the value in VALUE.  On return with STATUS zero, ARG is the routine's
c    estimate for the function's zero.
c
c    Input/output, integer STATUS, used to communicate between 
c    the user and the routine.  The user only sets STATUS to zero on the first 
c    call, to indicate that this is a startup call.  The routine returns STATUS
c    positive to request that the function be evaluated at ARG, or returns
c    STATUS as 0, to indicate that the iteration is complete and that
c    ARG is the estimated zero
c
c    Input, double precision VALUE, the function value at ARG, as requested
c    by the routine on the previous call.
c
      implicit none

      double precision a
      double precision arg
      double precision b
      double precision c
      save c
      double precision d
      save d
      double precision e
      save e
      double precision fa
      save fa
      double precision fb
      save fb
      double precision fc
      save fc
      double precision m
      double precision machep
      save machep
      double precision p
      double precision q
      double precision r
      double precision r8_epsilon
      double precision s
      double precision sa
      save sa
      double precision sb
      save sb
      integer status
      double precision t
      double precision tol
      double precision value
c
c  Input STATUS = 0.
c  Initialize, request F(A).
c
      if ( status .eq. 0 ) then

        machep = r8_epsilon ( a )

        sa = a
        sb = b
        e = sb - sa
        d = e

        status = 1
        arg = a
        return
c
c  Input STATUS = 1.
c  Receive F(A), request F(B).
c
      else if ( status .eq. 1 ) then

        fa = value

        status = 2
        arg = sb
        return
c
c  Input STATUS = 2
c  Receive F(B).
c
      else if ( status .eq. 2 ) then

        fb = value

        if ( 0.0D+00 .lt. fa * fb ) then
          status = -1
          return
        end if

        c = sa
        fc = fa

      else

        fb = value

        if ( ( 0.0D+00 .lt. fb .and. 0.0D+00 .lt. fc ) .or. 
     &       ( fb .le. 0.0D+00 .and. fc .le. 0.0D+00 ) ) then
          c = sa
          fc = fa
          e = sb - sa
          d = e
        end if

      end if
c
c  Compute the next point at which a function value is requested.
c
      if ( abs ( fc ) .lt. abs ( fb ) ) then

        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa

      end if

      tol = 2.0D+00 * machep * abs ( sb ) + t
      m = 0.5D+00 * ( c - sb )

      if ( abs ( m ) .le. tol .or. fb .eq. 0.0D+00 ) then
        status = 0
        arg = sb
        return
      end if

      if ( abs ( e ) .lt. tol .or. abs ( fa ) .le. abs ( fb ) ) then

        e = m
        d = e

      else

        s = fb / fa

        if ( sa .eq. c ) then

          p = 2.0D+00 * m * s
          q = 1.0D+00 - s

        else

          q = fa / fc
          r = fb / fc
          p = s * ( 2.0D+00 * m * q * ( q - r ) 
     &      - ( sb - sa ) * ( r - 1.0D+00 ) )
          q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

        end if

        if ( 0.0D+00 .lt. p ) then
          q = - q
        else
          p = - p
        end if

        s = e
        e = d

        if ( 2.0D+00 * p .lt. 3.0D+00 * m * q - abs ( tol * q ) .and. 
     &    p .lt. abs ( 0.5D+00 * s * q ) ) then
          d = p / q
        else
          e = m
          d = e
        end if

      end if

      sa = sb
      fa = fb

      if ( tol .lt. abs ( d ) ) then
        sb = sb + d
      else if ( 0.0D+00 .lt. m ) then
        sb = sb + tol
      else
        sb = sb - tol
      end if

      arg = sb
      status = status + 1

      return
      end
