      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8_epsilon

      r8_epsilon = 2.220446049250313D-016

      return
      end
      function r8_factorial ( n )

c*********************************************************************72
c
cc R8_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c
c    Output, double precision R8_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer n
      double precision r8_factorial

      r8_factorial = 1.0D+00

      do i = 1, n
        r8_factorial = r8_factorial * dble ( i )
      end do

      return
      end
      function r8_factorial2 ( n )

c*********************************************************************72
c
cc R8_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Example:
c
c     N   Value
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the double factorial 
c    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
c
c    Output, double precision R8_FACTORIAL2, the value.
c
      implicit none

      integer n
      double precision r8_factorial2
      double precision r8_n

      if ( n .lt. 1 ) then
        r8_factorial2 = 1.0D+00
        return
      end if

      r8_n = dble ( n )
      r8_factorial2 = 1.0D+00

10    continue

      if ( 1.0D+00 .lt. r8_n ) then
        r8_factorial2 = r8_factorial2 * r8_n
        r8_n = r8_n - 2.0D+00
        go to 10
      end if

      return
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates Gamma(X) for a real argument.
c
c  Discussion:
c
c    This routine calculates the gamma function for a real argument X.
c    Computation is based on an algorithm outlined in reference 1.
c    The program uses rational functions that approximate the gamma
c    function to at least 20 significant decimal digits.  Coefficients
c    for the approximation over the interval (1,2) are unpublished.
c    Those for the approximation for 12 <= X are from reference 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    An Overview of Software Development for Special Functions,
c    in Numerical Analysis Dundee, 1975,
c    edited by GA Watson,
c    Lecture Notes in Mathematics 506,
c    Springer, 1976.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision eps
      double precision fact
      integer i
      integer n
      double precision p(8)
      logical parity
      double precision r8_pi
      double precision q(8)
      double precision r8_gamma
      double precision res
      double precision sqrtpi
      double precision sum
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
c
c  Mathematical constants
c
      data sqrtpi /0.9189385332046727417803297D+00/
      data r8_pi /3.1415926535897932384626434D+00/
c
c  Machine dependent parameters
c
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps /2.22D-16/
      data xinf /1.79D+308/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (1,2).
c
      data p/
     & -1.71618513886549492533811d+00,
     &  2.47656508055759199108314d+01,
     & -3.79804256470945635097577d+02,
     &  6.29331155312818442661052d+02,
     &  8.66966202790413211295064d+02,
     & -3.14512729688483675254357d+04,
     & -3.61444134186911729807069d+04,
     &  6.64561438202405440627855d+04/

      data q/
     & -3.08402300119738975254353d+01,
     &  3.15350626979604161529144d+02,
     & -1.01515636749021914166146d+03,
     & -3.10777167157231109440444d+03,
     &  2.25381184209801510330112d+04,
     &  4.75584627752788110767815d+03,
     & -1.34659959864969306392456d+05,
     & -1.15132259675553483497211d+05/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     & -1.910444077728D-03,
     &  8.4171387781295D-04,
     & -5.952379913043012D-04,
     &  7.93650793500350248D-04,
     & -2.777777777777681622553D-03,
     &  8.333333333333333331554247D-02,
     &  5.7083835261D-03/

      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
c
c  Argument is negative.
c
      if ( y .le. 0.0D+00 ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. 0.0D+00 ) then

          if ( y1 .ne. aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
            parity = .true.
          end if

          fact = - r8_pi / sin ( r8_pi * res )
          y = y + 1.0D+00

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Argument is positive.
c
      if ( y .lt. eps ) then
c
c  Argument < EPS.
c
        if ( xminin .le. y ) then
          res = 1.0D+00 / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. 12.0D+00 ) then

        y1 = y
c
c  0.0 < argument < 1.0.
c
        if ( y .lt. 1.0D+00 ) then

          z = y
          y = y + 1.0D+00
c
c  1.0 < argument < 12.0.
c  Reduce argument if necessary.
c
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - 1.0D+00

        end if
c
c  Evaluate approximation for 1.0 < argument < 2.0.
c
        xnum = 0.0D+00
        xden = 1.0D+00
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + 1.0D+00
c
c  Adjust result for case  0.0 < argument < 1.0.
c
        if ( y1 .lt. y ) then

          res = res / y1
c
c  Adjust result for case 2.0 < argument < 12.0.
c
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + 1.0D+00
          end do

        end if

      else
c
c  Evaluate for 12.0 <= argument.
c
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - 0.5D+00 ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Final adjustments and return.
c
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. 1.0D+00 ) then
        res = fact / res
      end if

      r8_gamma = res

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function pythag ( a, b )

c*********************************************************************72
c
cc PYTHAG computes SQRT ( A * A + B * B ) carefully.
c
c  Discussion:
c
c    The formula
c
c      pythag = sqrt ( a * a + b * b )
c
c    is reasonably accurate, but can fail if, for example, A*A is larger
c    than the machine overflow.  The formula can lose most of its accuracy
c    if the sum of the squares is very large or very small.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 October 2009
c
c  Author:
c
c    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
c    Klema, Moler.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    James Wilkinson, Christian Reinsch,
c    Handbook for Automatic Computation,
c    Volume II, Linear Algebra, Part 2,
c    Springer, 1971,
c    ISBN: 0387054146,
c    LC: QA251.W67.
c
c    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
c    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
c    Matrix Eigensystem Routines, EISPACK Guide,
c    Lecture Notes in Computer Science, Volume 6,
c    Springer Verlag, 1976,
c    ISBN13: 978-3540075462,
c    LC: QA193.M37.
c
c  Parameters:
c
c    Input, double precision A, B, the two legs of a right triangle.
c
c    Output, double precision PYTHAG, the length of the hypotenuse.
c
      implicit none

      double precision a
      double precision b
      double precision p
      double precision pythag
      double precision r
      double precision s
      double precision t
      double precision u

      p = max ( abs ( a ), abs ( b ) )

      if ( p .ne. 0.0D+00 ) then

        r = ( min ( abs ( a ), abs ( b ) ) / p )**2

10      continue

          t = 4.0D+00 + r

          if ( t .eq. 4.0D+00 ) then
            go to 20
          end if

          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u ) * ( s / u ) * r

        go to 10

20      continue

      end if

      pythag = p

      return
      end
      subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

c*********************************************************************72
c
cc R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
c
c  Discussion:
c
c    A minor bug was corrected.  The HW variable, used in several places as
c    the "old" value of a quantity being iteratively improved, was not
c    being initialized.  JVB, 11 February 2008.
c
c    The original version of this program allowed the input arguments to
c    be modified, although they were restored to their input values before exit.
c    This is unacceptable if the input arguments are allowed to be constants.
c    The code has been modified so that the input arguments are never modified.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
c    This FORTRAN77 version by John Burkardt.
c
c    The original FORTRAN77 version of this routine is copyrighted by
c    Shanjie Zhang and Jianming Jin.  However, they give permission to
c    incorporate this routine into a user program provided that the copyright
c    is acknowledged.
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45
c
c  Parameters:
c
c    Input, double precision A_INPUT, B_INPUT, C_INPUT, X_INPUT, 
c    the arguments of the function.  The user is allowed to pass these
c    values as constants or variables.
c    C_INPUT must not be equal to a nonpositive integer.
c    X_INPUT .lt. 1.
c
c    Output, double precision HF, the value of the function.
c
      implicit none

      double precision a
      double precision a_input
      double precision a0
      double precision aa
      double precision b
      double precision b_input
      double precision bb
      double precision c
      double precision c_input
      double precision c0
      double precision c1
      double precision el
      parameter ( el = 0.5772156649015329D+00 )
      double precision eps
      double precision f0
      double precision f1
      double precision g0
      double precision g1
      double precision g2
      double precision g3
      double precision ga
      double precision gabc
      double precision gam
      double precision gb
      double precision gbm
      double precision gc
      double precision gca
      double precision gcab
      double precision gcb
      double precision gm
      double precision hf
      double precision hw
      integer j
      integer k
      logical l0
      logical l1
      logical l2
      logical l3
      logical l4
      logical l5
      integer m
      integer nm
      double precision pa
      double precision pb
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 )
      double precision r
      double precision r0
      double precision r1
      double precision r8_gamma
      double precision r8_psi
      double precision rm
      double precision rp
      double precision sm
      double precision sp
      double precision sp0
      double precision x
      double precision x_input
      double precision x1
c
c  Immediately copy the input argumentsc
c
      a = a_input
      b = b_input
      c = c_input
      x = x_input

      l0 = ( c .eq. aint ( c ) ) .and. ( c .lt. 0.0D+00 )
      l1 = ( 1.0D+00 - x .lt. 1.0D-15 ) .and. ( c - a - b .le. 0.0D+00 )
      l2 = ( a .eq. aint ( a ) ) .and. ( a .lt. 0.0D+00 )
      l3 = ( b .eq. aint ( b ) ) .and. ( b .lt. 0.0D+00 )
      l4 = ( c - a .eq. aint ( c - a ) ) .and. ( c - a .le. 0.0D+00 )
      l5 = ( c - b .eq. aint ( c - b ) ) .and. ( c - b .le. 0.0D+00 )

      if ( l0 .or. l1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
        write ( *, '(a)' ) '  The hypergeometric series is divergent.'
        return
      end if

      if ( 0.95D+00 .lt. x ) then
        eps = 1.0D-08
      else
        eps = 1.0D-15
      end if

      if ( x .eq. 0.0D+00 .or. a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then

        hf = 1.0D+00
        return

      else if ( 1.0D+00 - x .eq. eps .and. 0.0D+00 .lt. c - a - b ) then

        gc = r8_gamma ( c )
        gcab = r8_gamma ( c - a - b )
        gca = r8_gamma ( c - a )
        gcb = r8_gamma ( c - b )
        hf = gc * gcab / ( gca * gcb )
        return

      else if ( 1.0D+00 + x .le. eps .and. 
     &  abs ( c - a + b - 1.0D+00 ) .le. eps ) then

        g0 = sqrt ( r8_pi ) * 2.0D+00**( - a )
        g1 = r8_gamma ( c )
        g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
        g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
        hf = g0 * g1 / ( g2 * g3 )
        return

      else if ( l2 .or. l3 ) then

        if ( l2 ) then
          nm = int ( abs ( a ) )
        end if

        if ( l3 ) then
          nm = int ( abs ( b ) )
        end if

        hf = 1.0D+00
        r = 1.0D+00

        do k = 1, nm
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do

        return

      else if ( l4 .or. l5 ) then

        if ( l4 ) then
          nm = int ( abs ( c - a ) )
        end if

        if ( l5 ) then
          nm = int ( abs ( c - b ) )
        end if

        hf = 1.0D+00
        r  = 1.0D+00
        do k = 1, nm
          r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        hf = ( 1.0D+00 - x )**( c - a - b ) * hf
        return

      end if

      aa = a
      bb = b
      x1 = x

      if ( x .lt. 0.0D+00 ) then
        x = x / ( x - 1.0D+00 )
        if ( a .lt. c .and. b .lt. a .and. 0.0D+00 .lt. b ) then
          a = bb
          b = aa
        end if
        b = c - b
      end if

      if ( 0.75D+00 .le. x ) then

        gm = 0.0D+00

        if ( abs ( c - a - b - aint ( c - a - b ) ) .lt. 1.0D-15 ) then

          m = int ( c - a - b )
          ga = r8_gamma ( a )
          gb = r8_gamma ( b )
          gc = r8_gamma ( c )
          gam = r8_gamma ( a + m )
          gbm = r8_gamma ( b + m )

          pa = r8_psi ( a )
          pb = r8_psi ( b )

          if ( m /= 0 ) then
            gm = 1.0D+00
          end if

          do j = 1, abs ( m ) - 1
            gm = gm * j
          end do

          rm = 1.0D+00
          do j = 1, abs ( m )
            rm = rm * j
          end do

          f0 = 1.0D+00
          r0 = 1.0D+00
          r1 = 1.0D+00
          sp0 = 0.0D+00
          sp = 0.0D+00

          if ( 0 .le. m ) then

            c0 = gm * gc / ( gam * gbm )
            c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

            do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) 
     &          + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / dble ( k )
            end do

            f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1

            do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) 
     &          + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                sm = sm + ( 1.0D+00 - a ) 
     &            / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) 
     &            + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do

              rp = pa + pb + 2.0D+00 * el + sp + sm 
     &          + log ( 1.0D+00 - x )

              r1 = r1 * ( a + m + k - 1.0D+00 ) 
     &          * ( b + m + k - 1.0D+00 ) 
     &          / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) .lt. abs ( f1 ) * eps ) then
                exit
              end if

              hw = f1

            end do

            hf = f0 * c0 + f1 * c1

          else if ( m .lt. 0 ) then

            m = - m
            c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
            c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

            do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) 
     &          * ( b - m + k - 1.0D+00 ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / dble ( k )
            end do

            f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1

            do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) 
     &          / ( k * ( a + k - 1.0D+00 ) ) 
     &          + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                sm = sm + 1.0D+00 / dble ( j + k )
              end do

              rp = pa + pb + 2.0D+00 * el + sp - sm 
     &          + log ( 1.0D+00 - x )

              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &          / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) .lt. abs ( f1 ) * eps ) then
                exit
              end if

              hw = f1

            end do

            hf = f0 * c0 + f1 * c1

          end if

        else

          ga = r8_gamma ( a )
          gb = r8_gamma ( b )
          gc = r8_gamma ( c )
          gca = r8_gamma ( c - a )
          gcb = r8_gamma ( c - b )
          gcab = r8_gamma ( c - a - b )
          gabc = r8_gamma ( a + b - c )
          c0 = gc * gcab / ( gca * gcb )
          c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
          hf = 0.0D+00
          hw = hf
          r0 = c0
          r1 = c1

          do k = 1, 250

            r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &        / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

            r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 )
     &        / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

            hf = hf + r0 + r1

            if ( abs ( hf - hw ) .lt. abs ( hf ) * eps ) then
              exit
            end if

            hw = hf

          end do

          hf = hf + c0 + c1

        end if

      else

        a0 = 1.0D+00

        if ( a .lt. c .and. c .lt. 2.0D+00 * a .and. 
     &       b .lt. c .and. c .lt. 2.0D+00 * b ) then

          a0 = ( 1.0D+00 - x )**( c - a - b )
          a = c - a
          b = c - b

        end if

        hf = 1.0D+00
        hw = hf
        r = 1.0D+00

        do k = 1, 250

          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x

          hf = hf + r

          if ( abs ( hf - hw ) .le. abs ( hf ) * eps ) then
            exit
          end if

          hw = hf

        end do

        hf = a0 * hf

      end if

      if ( x1 .lt. 0.0D+00 ) then
        x = x1
        c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
        hf = c0 * hf
      end if

      a = aa
      b = bb

      if ( 120 .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
        write ( *, '(a)' ) '  A large number of iterations were needed.'
        write ( *, '(a)' ) 
     &    '  The accuracy of the results should be checked.'
      end if

      return
      end
      function r8_psi ( xx )

c*********************************************************************72
c
cc R8_PSI evaluates the function Psi(X).
c
c  Discussion:
c
c    This routine evaluates the logarithmic derivative of the
c    GAMMA function,
c
c      PSI(X) = d/dX (GAMMA(X)) / GAMMA(X) 
c             = d/dX LN ( GAMMA(X) )
c
c    for real X, where either
c
c      -XMAX1 < X < -XMIN  and X is not a negative integer), 
c
c    or
c
c      XMIN < X.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, Anthony Strecok, Henry Thacher,
c    Chebyshev Approximations for the Psi Function,
c    Mathematics of Computation,
c    Volume 27, Number 121, January 1973, pages 123-127.
c
c  Parameters:
c
c    Input, double precision XX, the argument of the function.
c
c    Output, double precision R8_PSI, the value of the function.
c
      implicit none

      double precision aug
      double precision den
      integer i
      integer n
      integer nq
      double precision p1(9)
      double precision p2(7)
      double precision piov4
      double precision q1(8)
      double precision q2(6)
      double precision r8_psi
      double precision sgn
      double precision xlarge
      double precision upper
      double precision w
      double precision x
      double precision xinf
      double precision xmax1
      double precision xmin1
      double precision xsmall
      double precision x01
      double precision x01d
      double precision x02
      double precision xx
      double precision z
c
c  Mathematical constants.  PIOV4 = pi / 4
c
      data piov4 /7.8539816339744830962d-01/
c
c  Machine-dependent constants
c
      data xinf /1.70d+38/
      data xmin1 /5.89d-39/
      data xmax1 /3.60d+16/
      data xsmall /2.05d-09/
      data xlarge /2.04d+15/
c
c  Zero of psi(x)
c
      data x01 /187.0d0/
      data x01d /128.0d0/
      data x02 /6.9464496836234126266d-04/
c
c  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
c
      data p1/4.5104681245762934160d-03,5.4932855833000385356d+00,
     &        3.7646693175929276856d+02,7.9525490849151998065d+03,
     &        7.1451595818951933210d+04,3.0655976301987365674d+05,
     &        6.3606997788964458797d+05,5.8041312783537569993d+05,
     &        1.6585695029761022321d+05/
      data q1/9.6141654774222358525d+01,2.6287715790581193330d+03,
     &        2.9862497022250277920d+04,1.6206566091533671639d+05,
     &        4.3487880712768329037d+05,5.4256384537269993733d+05,
     &        2.4242185002017985252d+05,6.4155223783576225996d-08/
c
c  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
c  for 3.0 < x.
c
      data p2/-2.7103228277757834192d+00,-1.5166271776896121383d+01,
     &        -1.9784554148719218667d+01,-8.8100958828312219821d+00,
     &        -1.4479614616899842986d+00,-7.3689600332394549911d-02,
     &        -6.5135387732718171306d-21/
      data q2/ 4.4992760373789365846d+01, 2.0240955312679931159d+02,
     &         2.4736979003315290057d+02, 1.0742543875702278326d+02,
     &         1.7463965060678569906d+01, 8.8427520398873480342d-01/

      x = xx
      w = abs ( x )
      aug = 0.0D+00
c
c  Check for valid arguments, then branch to appropriate algorithm.
c
      if ( - x .ge. xmax1 .or. w .lt. xmin1 ) then
        r8_psi = xinf
        if ( 0.0D+00 .lt. x ) then
          r8_psi = -xinf
        end if
        return
      end if

      if ( x .ge. 0.5D+00 ) then
        go to 200
c
c  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
c  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
c
      else if ( w .le. xsmall ) then
        aug = - 1.0D+00 / x
        go to 150
      end if
c
c  Argument reduction for cotangent.
c
  100 continue

      if ( x .lt. 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - aint ( w )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - dble ( nq ) * 0.25D+00 )
c
c  W is now related to the fractional part of 4.0 * X.
c  Adjust argument to correspond to values in the first
c  quadrant and determine the sign.
c
      n = nq / 2

      if ( n + n .ne. nq ) then
        w = 1.0D+00 - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) .ne. 0 ) then
        sgn = - sgn
      end if
c
c  Determine the final value for  -pi * cotan(pi*x).
c
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) .eq. 0 ) then
c
c  Check for singularity.
c
        if ( z .eq. 0.0D+00 ) then
          r8_psi = xinf
          if ( 0.0D+00 .lt. x ) then
            r8_psi = -xinf
          end if
          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else
        aug = sgn * ( 4.0D+00 * tan ( z ) )
      end if

  150 continue

      x = 1.0D+00 - x

  200 continue
c
c  0.5 <= X <= 3.0.
c
      if ( x .le. 3.0D+00 ) then

        den = x
        upper = p1(1) * x
        do i = 1, 7
          den = ( den + q1(i) ) * x
          upper = ( upper + p1(i+1) ) * x
        end do
        den = ( upper + p1(9) ) / ( den + q1(8) )
        x = ( x - x01 / x01d ) - x02
        r8_psi = den * x + aug
        return

      end if
c
c  3.0 < X.
c
      if ( x .lt. xlarge ) then
        w = 1.0D+00 / ( x * x )
        den = w
        upper = p2(1) * w
        do i = 1, 5
          den = ( den + q2(i) ) * w
          upper = ( upper + p2(i+1) ) * w
        end do
        aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
      end if

      r8_psi = aug + log ( x )

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      function r8vec_product ( n, v1 )

c*********************************************************************72
c
cc R8VEC_PRODUCT multiplies the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_PRODUCT, the product of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_product
      double precision v1(n)
      double precision value

      value = 1.0D+00
      do i = 1, n
        value = value * v1(i)
      end do

      r8vec_product = value

      return
      end
      subroutine r8vec_reverse ( n, a )

c*********************************************************************72
c
cc R8VEC_REVERSE reverses the elements of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, calling R8VEC_REVERSE is equivalent to
c
c      A(1:N) = A(N:1:-1)
c
c  Example:
c
c    Input:
c
c      N = 5,
c      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
c
c    Output:
c
c      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, double precision A(N), the array to be reversed.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision t

      do i = 1, n / 2
        t        = a(i)
        a(i)     = a(n+1-i)
        a(n+1-i) = t
      end do

      return
      end
