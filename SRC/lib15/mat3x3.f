      subroutine mat3x3( a, b )
      implicit none
c based on routine GAUSSJ of "Numerical Recipes"
c
c calling arguments
      real*8 a( 3, 3 ), b( 3 )
c
c local arrays
      integer indxc( 3 ), indxr( 3 ), ipiv( 3 )
c
c local variables
      integer i, icol, irow, j, k, l, ll
      real*8 big, dum, pivinv
c
      do j = 1, 3
        ipiv( j ) = 0
      end do
c loop over columns
      do i = 1, 3
c search for the pivot element
        big = 0.
        do j = 1, 3
          if( ipiv( j ) .ne. 1 )then
            do k = 1, 3
              if( ipiv( k ) .eq. 0 )then
                if( abs( a( j, k ) ) .ge. big )then
                  big = abs( a( j, k ) )
                  irow = j
                  icol = k
                end if
              else if( ipiv( k ) .gt. 1 )then
                call crash( 'MAT3X3', 'singular matrix 1' )
              end if
            end do
          end if
        end do
        ipiv( icol ) = ipiv( icol ) + 1
c interchange rows if needed
        if( irow .ne. icol )then
          do l = 1, 3
            dum = a( irow, l )
            a( irow, l ) = a( icol, l )
            a( icol, l ) = dum
          end do
          do l = 1, 1
            dum = b( irow )
            b( irow ) = b( icol )
            b( icol ) = dum
          end do
        end if
c divide by the selected pivot element
        indxr( i ) = irow
        indxc( i ) = icol
        if( a( icol, icol )
     +              .eq. 0. )call crash( 'MAT3X3', 'singular matrix 2' )
        pivinv = 1. / a( icol, icol )
        a( icol, icol ) = 1.
        do l = 1, 3
          a( icol, l ) = a( icol, l ) * pivinv
        end do
        b( icol ) = b( icol ) * pivinv
c reduce the rows
        do ll = 1, 3
          if( ll .ne. icol )then
            dum = a( ll, icol )
            a( ll, icol ) = 0.
            do l = 1, 3
              a( ll, l ) = a( ll, l ) - a( icol, l ) * dum
            end do
            b( ll ) = b( ll ) - b( icol ) * dum
          end if
        end do
c end loop over columns
      end do
c unscramble the solution
      do l = 3, 1, -1
        if( indxr( l ) .ne. indxc( l ) )then
          do k = 1, 3
            dum = a( k, indxr( l ) )
            a( k, indxr( l ) ) = a( k, indxc( l ) )
            a( k, indxc( l ) ) = dum
          end do
        end if
      end do
      return
      end
