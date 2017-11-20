      subroutine imtql2 ( n, d, e, z, ierr )

c*********************************************************************72
c
cc IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
c
c  Discussion:
c
c    This subroutine finds the eigenvalues and eigenvectors
c    of a symmetric tridiagonal matrix by the implicit QL method.
c    The eigenvectors of a full symmetric matrix can also
c    be found if TRED2 has been used to reduce this
c    full matrix to tridiagonal form.
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
c    FORTRAN90 version by John Burkardt.
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
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision D(N).  On input, the diagonal elements of
c    the input matrix.  On output, the eigenvalues in ascending order.  If an
c    error exit is made, the eigenvalues are correct but
c    unordered for indices 1,2,...,IERR-1.
c
c    Input/output, double precision E(N).  On input, the subdiagonal elements
c    of the input matrix in E(2:N).  E(1) is arbitrary.  On output, E is
c    overwritten.
c
c    Input/output, double precision Z(N,N).  On input, the transformation
c    matrix produced in the reduction by TRED2, if performed.  If the
c    eigenvectors of the tridiagonal matrix are desired, Z must contain the
c    identity matrix.  On output, Z contains orthonormal eigenvectors of the
c    symmetric tridiagonal (or full) matrix.  If an error exit is made, Z
c    contains the eigenvectors associated with the stored eigenvalues.
c
c    Output, integer IERR, error flag.
c    0, for normal return,
c    J, if the J-th eigenvalue has not been determined after 30 iterations.
c
      implicit none

      integer n

      double precision b
      double precision c
      double precision d(n)
      double precision e(n)
      double precision f
      double precision g
      integer i
      integer ierr
      integer ii
      integer j
      integer k
      integer l
      integer m
      integer mml
      double precision p
      double precision pythag
      double precision r
      double precision s
      double precision t(n)
      double precision tst1
      double precision tst2
      double precision z(n,n)

      ierr = 0

      if ( n .eq. 1 ) then
        return
      end if

      do i = 2, n
        e(i-1) = e(i)
      end do
      e(n) = 0.0D+00

      do l = 1, n

        j = 0
c
c  Look for a small sub-diagonal element.
c
105     continue

        do m = l, n

          if ( m .eq. n ) then
            go to 106
          end if

          tst1 = abs ( d(m) ) + abs ( d(m+1) )
          tst2 = tst1 + abs ( e(m) )

          if ( tst2 .eq. tst1 ) then
            go to 106
          end if

        end do

106     continue

        p = d(l)

        if ( m .eq. l ) then
          go to 107
        end if

        if ( 30 .le. j ) then
          ierr = l
          return
        end if

        j = j + 1
c
c  Form shift.
c
        g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
        r = pythag ( g, 1.0D+00 )
        g = d(m) - p + e(l) / ( g + sign ( r, g ) )
        s = 1.0D+00
        c = 1.0D+00
        p = 0.0D+00
        mml = m - l

        do ii = 1, mml

          i = m - ii
          f = s * e(i)
          b = c * e(i)
          r = pythag ( f, g )
          e(i+1) = r
c
c  Recover from underflow.
c
          if ( r .eq. 0.0D+00 ) then
            d(i+1) = d(i+1) - p
            e(m) = 0.0D+00
            go to 105
          end if

          s = f / r
          c = g / r
          g = d(i+1) - p
          r = ( d(i) - g ) * s + 2.0D+00 * c * b
          p = s * r
          d(i+1) = g + p
          g = c * r - b
c
c  Form vector.
c
          do k = 1, n
            f = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * f
            z(k,i) = c * z(k,i) - s * f
          end do

        end do

        d(l) = d(l) - p
        e(l) = g
        e(m) = 0.0D+00
        go to 105

107     continue

      end do
c
c  Order eigenvalues and eigenvectors.
c
      do ii = 2, n

        i = ii - 1
        k = i
        p = d(i)

        do j = ii, n
          if ( d(j) .lt. p ) then
            k = j
            p = d(j)
          end if
        end do

        if ( k .ne. i ) then

          d(k) = d(i)
          d(i) = p

          do j = 1, n
            t(j)   = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = t(j)
          end do

        end if

      end do

      return
      end
      subroutine imtqlx ( n, d, e, z )

c*********************************************************************72
c
cc IMTQLX diagonalizes a symmetric tridiagonal matrix.
c
c  Discussion:
c
c    This routine is a slightly modified version of the EISPACK routine to
c    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
c
c    The authors thank the authors of EISPACK for permission to use this
c    routine.
c
c    It has been modified to produce the product Q' * Z, where Z is an input
c    vector and Q is the orthogonal matrix diagonalizing the input matrix.
c    The changes consist (essentially) of applying the orthogonal
c    transformations directly to Z as they are generated.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 April 2011
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
c    Roger Martin, James Wilkinson,
c    The Implicit QL Algorithm,
c    Numerische Mathematik,
c    Volume 12, Number 5, December 1968, pages 377-383.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision D(N), the diagonal entries of the matrix.
c    On output, the information in D has been overwritten.
c
c    Input/output, double precision E(N), the subdiagonal entries of the
c    matrix, in entries E(1) through E(N-1).  On output, the information in
c    E has been overwritten.
c
c    Input/output, double precision Z(N).  On input, a vector.  On output,
c    the value of Q' * Z, where Q is the matrix that diagonalizes the
c    input symmetric tridiagonal matrix.
c
      implicit none

      integer n

      double precision b
      double precision c
      double precision d(n)
      double precision e(n)
      double precision f
      double precision g
      integer i
      integer ii
      integer itn
      parameter ( itn = 30 )
      integer j
      integer k
      integer l
      integer m
      integer mml
      double precision p
      double precision prec
      double precision r
      double precision r8_epsilon
      double precision s
      double precision test
      double precision z(n)

      prec = r8_epsilon ( prec )

      if ( n .eq. 1 ) then
        return
      end if

      e(n) = 0.0D+00

      do l = 1, n

        j = 0

10      continue

          do m = l, n

            if ( m == n ) then
              go to 20
            end if

            test = prec * ( abs ( d(m) ) + abs ( d(m+1) ) )

            if ( abs ( e(m) ) .le. test ) then
              go to 20
            end if

          end do

20        continue

          p = d(l)

          if ( m .eq. l ) then
            go to 30
          end if

          if ( itn .le. j ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'IMTQLX - Fatal error!'
            write ( *, '(a)' ) '  Iteration limit exceeded.'
            write ( *, '(a,i8)' ) '  J = ', j
            write ( *, '(a,i8)' ) '  L = ', l
            write ( *, '(a,i8)' ) '  M = ', m
            write ( *, '(a,i8)' ) '  N = ', n
            stop 1
          end if

          j = j + 1
          g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
          r =  sqrt ( g * g + 1.0D+00 )
          g = d(m) - p + e(l) / ( g + sign ( r, g ) )
          s = 1.0D+00
          c = 1.0D+00
          p = 0.0D+00
          mml = m - l

          do ii = 1, mml

            i = m - ii
            f = s * e(i)
            b = c * e(i)

            if ( abs ( g ) .le. abs ( f ) ) then
              c = g / f
              r =  sqrt ( c * c + 1.0D+00 )
              e(i+1) = f * r
              s = 1.0D+00 / r
              c = c * s
            else
              s = f / g
              r =  sqrt ( s * s + 1.0D+00 )
              e(i+1) = g * r
              c = 1.0D+00 / r
              s = s * c
            end if

            g = d(i+1) - p
            r = ( d(i) - g ) * s + 2.0D+00 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f

          end do

          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0D+00

        go to 10

30      continue

      end do
c
c  Sorting.
c
      do ii = 2, n

        i = ii - 1
        k = i
        p = d(i)

        do j = ii, n
          if ( d(j) .lt. p ) then
            k = j
            p = d(j)
          end if
        end do

        if ( k .ne. i ) then
          d(k) = d(i)
          d(i) = p
          p = z(i)
          z(i) = z(k)
          z(k) = p
        end if

      end do

      return
      end
