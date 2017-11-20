      subroutine jacobi_ek_compute ( n, alpha, beta, x, w )

c*********************************************************************72
c
cc JACOBI_EK_COMPUTE: Elhay-Kautsky method for Gauss-Jacobi quadrature rule.
c
c  Discussion:
c
c    The integral:
c
c      integral ( -1 <= x <= 1 ) (1-x)^alpha * (1+x)^beta * f(x) dx
c
c    The quadrature rule:
c
c      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 April 2011
c
c  Author:
c
c    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Sylvan Elhay, Jaroslav Kautsky,
c    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
c    Interpolatory Quadrature,
c    ACM Transactions on Mathematical Software,
c    Volume 13, Number 4, December 1987, pages 399-415.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision ALPHA, BETA, the exponents of (1-X) and
c    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
c    set ALPHA = BETA = 0.0.  -1.0 .lt. ALPHA and -1.0 .lt. BETA are required.
c
c    Output, double precision X(N), the abscissas.
c
c    Output, double precision W(N), the weights.
c
      implicit none

      integer n

      double precision alpha
      double precision abi
      double precision beta
      double precision bj(n)
      integer i
      double precision i_r8
      double precision r8_gamma
      double precision w(n)
      double precision x(n)
      double precision zemu
c
c  Define the zero-th moment.
c
      zemu = 2.0D+00**( alpha + beta + 1.0D+00 ) 
     &  * r8_gamma ( alpha + 1.0D+00 ) 
     &  * r8_gamma ( beta + 1.0D+00 ) 
     &  / r8_gamma ( 2.0D+00 + alpha + beta )
c
c  Define the Jacobi matrix.
c
      x(1) = ( beta - alpha ) / ( 2.0D+00 + alpha + beta )

      bj(1) = 4.0D+00 * ( 1.0D+00 + alpha ) * ( 1.0D+00 + beta ) 
     &  / ( ( 3.0D+00 + alpha + beta ) * ( 2.0D+00 + alpha + beta )**2 )

      do i = 2, n
        i_r8 = dble ( i )
        abi = 2.0D+00 * i_r8 + alpha + beta
        x(i) = ( beta + alpha ) * ( beta - alpha ) 
     &    / ( ( abi - 2.0D+00 ) * abi )
        bj(i) = 4.0D+00 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) 
     &    * ( i_r8 + alpha + beta ) 
     &    / ( ( abi - 1.0D+00 ) * ( abi + 1.0D+00 ) * abi * abi )
      end do

      do i = 1, n
        bj(i) = sqrt ( bj(i) )
      end do

      w(1) = sqrt ( zemu )
      do i = 2, n 
        w(i) = 0.0D+00
      end do
c
c  Diagonalize the Jacobi matrix.
c
      call imtqlx ( n, x, bj, w )

      do i = 1, n
        w(i) = w(i)**2
      end do

      return
      end
      subroutine jacobi_integral ( expon, alpha, beta, value )

c*********************************************************************72
c
cc JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
c
c  Discussion:
c
c    The integral:
c
c      integral ( -1 <= x <= +1 ) x^expon (1-x)^alpha (1+x)^beta dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer EXPON, the exponent.
c
c    Input, double precision ALPHA, the exponent of (1-X) in the weight factor.
c
c    Input, double precision BETA, the exponent of (1+X) in the weight factor.
c
c    Output, double precision VALUE, the value of the integral.
c
      implicit none

      double precision alpha
      double precision arg1
      double precision arg2
      double precision arg3
      double precision arg4
      double precision beta
      double precision c
      integer expon
      double precision r8_gamma
      double precision s
      double precision value
      double precision value1
      double precision value2

      c = dble ( expon )

      if ( mod ( expon, 2 ) .eq. 0 ) then
        s = +1.0D+00
      else
        s = -1.0D+00
      end if

      arg1 = - alpha
      arg2 =   1.0D+00 + c
      arg3 =   2.0D+00 + beta + c
      arg4 = - 1.0D+00

      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

      arg1 = - beta
      arg2 =   1.0D+00 + c
      arg3 =   2.0D+00 + alpha + c
      arg4 = - 1.0D+00

      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )

      value = r8_gamma ( 1.0D+00 + c ) * ( 
     &    s * r8_gamma ( 1.0D+00 + beta  ) * value1 
     &  / r8_gamma ( 2.0D+00 + beta  + c ) 
     &  +     r8_gamma ( 1.0D+00 + alpha ) * value2 
     &  / r8_gamma ( 2.0D+00 + alpha + c ) )

      return
      end
      subroutine jacobi_ss_compute ( n, alpha, beta, x, w )

c*********************************************************************72
c
cc JACOBI_SS_COMPUTE computes a Gauss-Jacobi quadrature rule.
c
c  Discussion:
c
c    The integral:
c
c      integral ( -1 .le. x .le. 1 ) (1-x)^alpha * (1+x)^beta * f(x) dx
c
c    The quadrature rule:
c
c      sum ( 1 .le. i .le. n ) w(i) * f ( x(i) )
c
c    Thanks to Xu Xiang of Fudan University for pointing out that
c    an earlier implementation of this routine was incorrect!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 May 2007
c
c  Author:
c
c    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision ALPHA, BETA, the exponents of (1-X) and
c    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
c    set ALPHA = BETA = 0.0.  -1.0 .lt. ALPHA and -1.0 .lt. BETA are required.
c
c    Output, double precision X(N), the abscissas.
c
c    Output, double precision W(N), the weights.
c
      implicit none

      integer n

      double precision alpha
      double precision an
      double precision b(n)
      double precision beta
      double precision bn
      double precision c(n)
      double precision cc
      double precision delta
      double precision dp2
      integer i
      double precision p1
      double precision r1
      double precision r2
      double precision r3
      double precision r8_gamma
      double precision r8vec_product
      double precision w(n)
      double precision x(n)
      double precision xval
c
c  Check ALPHA and BETA.
c
      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
        write ( *, '(a)' ) '  -1.0 .lt. ALPHA is required.'
        stop 1
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_SS_COMPUTE - Fatal error!'
        write ( *, '(a)' ) '  -1.0 .lt. BETA is required.'
        stop 1
      end if
c
c  Set the recursion coefficients.
c
      do i = 1, n

        if ( alpha + beta .eq. 0.0D+00 .or. 
     &    beta - alpha .eq. 0.0D+00 ) then

          b(i) = 0.0D+00

        else

          b(i) = ( alpha + beta ) * ( beta - alpha ) / 
     &          ( ( alpha + beta + dble ( 2 * i ) ) 
     &          * ( alpha + beta + dble ( 2 * i - 2 ) ) )

        end if

        if ( i .eq. 1 ) then

          c(i) = 0.0D+00

        else

          c(i) = 4.0D+00 * dble ( i - 1 ) 
     &          * ( alpha + dble ( i - 1 ) ) 
     &          * ( beta + dble ( i - 1 ) ) 
     &          * ( alpha + beta + dble ( i - 1 ) ) / 
     &          ( ( alpha + beta + dble ( 2 * i - 1 ) ) 
     &          * ( alpha + beta + dble ( 2 * i - 2 ) )**2 
     &          * ( alpha + beta + dble ( 2 * i - 3 ) ) )

        end if

      end do

      delta = r8_gamma ( alpha        + 1.0D+00 ) 
     &      * r8_gamma (         beta + 1.0D+00 ) 
     &      / r8_gamma ( alpha + beta + 2.0D+00 )

      cc = delta * 2.0D+00**( alpha + beta + 1.0D+00 ) 
     &  * r8vec_product ( n - 1, c(2) )

      do i = 1, n

        if ( i .eq. 1 ) then

          an = alpha / dble ( n )
          bn = beta / dble ( n )

          r1 = ( 1.0D+00 + alpha ) 
     &      * ( 2.78D+00 / ( 4.0D+00 + dble ( n * n ) ) 
     &      + 0.768D+00 * an / dble ( n ) )

          r2 = 1.0D+00 + 1.48D+00 * an + 0.96D+00 * bn 
     &      + 0.452D+00 * an * an + 0.83D+00 * an * bn

          xval = ( r2 - r1 ) / r2

        else if ( i .eq. 2 ) then

          r1 = ( 4.1D+00 + alpha ) / 
     &      ( ( 1.0D+00 + alpha ) * ( 1.0D+00 + 0.156D+00 * alpha ) )

          r2 = 1.0D+00 + 0.06D+00 * ( dble ( n ) - 8.0D+00 ) * 
     &      ( 1.0D+00 + 0.12D+00 * alpha ) / dble ( n )

          r3 = 1.0D+00 + 0.012D+00 * beta * 
     &      ( 1.0D+00 + 0.25D+00 * abs ( alpha ) ) / dble ( n )

          xval = xval - r1 * r2 * r3 * ( 1.0D+00 - xval )

        else if ( i .eq. 3 ) then

          r1 = ( 1.67D+00 + 0.28D+00 * alpha ) 
     &      / ( 1.0D+00 + 0.37D+00 * alpha )

          r2 = 1.0D+00 + 0.22D+00 * ( dble ( n ) - 8.0D+00 ) 
     &      / dble ( n )

          r3 = 1.0D+00 + 8.0D+00 * beta / 
     &      ( ( 6.28D+00 + beta ) * dble ( n * n ) )

          xval = xval - r1 * r2 * r3 * ( x(1) - xval )

        else if ( i .lt. n - 1 ) then

          xval = 3.0D+00 * x(i-1) - 3.0D+00 * x(i-2) + x(i-3)

        else if ( i .eq. n - 1 ) then

          r1 = ( 1.0D+00 + 0.235D+00 * beta ) 
     &      / ( 0.766D+00 + 0.119D+00 * beta )

          r2 = 1.0D+00 / ( 1.0D+00 + 0.639D+00 
     &      * ( dble ( n ) - 4.0D+00 ) 
     &      / ( 1.0D+00 + 0.71D+00 * ( dble ( n ) - 4.0D+00 ) ) )

          r3 = 1.0D+00 / ( 1.0D+00 + 20.0D+00 * alpha 
     &      / ( ( 7.5D+00 + alpha ) * dble ( n * n ) ) )

          xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

        else if ( i .eq. n ) then

          r1 = ( 1.0D+00 + 0.37D+00 * beta ) 
     &      / ( 1.67D+00 + 0.28D+00 * beta )

          r2 = 1.0D+00 / 
     &      ( 1.0D+00 + 0.22D+00 * ( dble ( n ) - 8.0D+00 ) 
     &      / dble ( n ) )

          r3 = 1.0D+00 / ( 1.0D+00 + 8.0D+00 * alpha / 
     &      ( ( 6.28D+00 + alpha ) * dble ( n * n ) ) )

          xval = xval + r1 * r2 * r3 * ( xval - x(i-2) )

        end if

        call jacobi_ss_root ( xval, n, alpha, beta, dp2, p1, b, c )

        x(i) = xval
        w(i) = cc / ( dp2 * p1 )

      end do
c
c  Reverse the data.
c
      call r8vec_reverse ( n, x )
      call r8vec_reverse ( n, w )

      return
      end
      subroutine jacobi_ss_recur ( p2, dp2, p1, x, n, alpha, beta, 
     &  b, c )

c*********************************************************************72
c
cc JACOBI_SS_RECUR finds the value and derivative of a Jacobi polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 1998
c
c  Author:
c
c    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c  Parameters:
c
c    Output, double precision P2, the value of J(N)(X).
c
c    Output, double precision DP2, the value of J'(N)(X).
c
c    Output, double precision P1, the value of J(N-1)(X).
c
c    Input, double precision X, the point at which polynomials are evaluated.
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision ALPHA, BETA, the exponents of (1-X) and
c    (1+X) in the quadrature rule.
c
c    Input, double precision B(N), C(N), the recursion coefficients.
c
      implicit none

      integer n

      double precision alpha
      double precision b(n)
      double precision beta
      double precision c(n)
      double precision dp0
      double precision dp1
      double precision dp2
      integer i
      double precision p0
      double precision p1
      double precision p2
      double precision x

      p1 = 1.0D+00
      dp1 = 0.0D+00

      p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0D+00 )
      dp2 = 1.0D+00

      do i = 2, n

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2 = ( x - b(i) ) * p1 - c(i) * p0
        dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

      end do

      return
      end
      subroutine jacobi_ss_root ( x, n, alpha, beta, dp2, p1, b, c )

c*********************************************************************72
c
cc JACOBI_SS_ROOT improves an approximate root of a Jacobi polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2000
c
c  Author:
c
c    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Arthur Stroud, Don Secrest,
c    Gaussian Quadrature Formulas,
c    Prentice Hall, 1966,
c    LC: QA299.4G3S7.
c
c  Parameters:
c
c    Input/output, double precision X, the approximate root, which
c    should be improved on output.
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision ALPHA, BETA, the exponents of (1-X) and
c    (1+X) in the quadrature rule.
c
c    Output, double precision DP2, the value of J'(N)(X).
c
c    Output, double precision P1, the value of J(N-1)(X).
c
c    Input, double precision B(N), C(N), the recursion coefficients.
c
      implicit none

      integer n

      double precision alpha
      double precision b(n)
      double precision beta
      double precision c(n)
      double precision d
      double precision dp2
      double precision eps
      double precision p1
      double precision p2
      double precision r8_epsilon
      integer step
      integer step_max
      parameter ( step_max = 10 )
      double precision x

      eps = r8_epsilon ( )

      do step = 1, step_max

        call jacobi_ss_recur ( p2, dp2, p1, x, n, alpha, beta, b, c )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + 1.0D+00 ) ) then
          return
        end if

      end do

      return
      end
