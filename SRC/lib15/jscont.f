      subroutine jscont( a, idim, jdim, c, nc, tr )
c  Copyright (c) 2014, Jerry Sellwood
      implicit none
c primitive routine to draw contours of the values specified in the
c    input array a at the nc levels specified in the input array c
c----------------------------------------------------------------------
c arguments:
c a (input, real array): data array.
c idim (input, integer): first dimension of a.
c jdim (input, integer): second dimension of a.
c c (input, real array): array of contour levels.
c nc (input, integer): number of contour levels (less than or
c       equal to dimension of c).
c tr (external function): transforms grid values to x,y screen
c       coordinates
c----------------------------------------------------------------------
c
c calling arguments
      external tr
      integer idim, jdim, nc
      real a( idim, jdim ), c( nc )
c
c local arrays
      integer id( 6 )
      real v( 5 ), x( 2 ), y( 2 )
c
c local variables
      integer i, ic, icorn, j, npt
      real ctr, delta, xx, yy
c
      data id / 0, -1, -1, 0, 0, -1 /
c
      if( nc .lt. 1 )return
c work over array
      do j = 2, jdim
        do i = 2, idim
c copy cell corner values
          v( 1 ) = a( i - 1, j )
          v( 2 ) = a( i - 1, j - 1 )
          v( 3 ) = a( i, j - 1 )
          v( 4 ) = a( i, j )
          v( 5 ) = v( 1 )
c scan over contours
          do ic = 1, nc
            ctr = c( ic )
            npt = 0
            do icorn = 1, 4
              if(
     +      ( v( icorn ) .ge. ctr .and. v( icorn + 1 ) .lt. ctr ) .or.
     +      ( v( icorn ) .lt. ctr .and. v( icorn + 1 ) .ge. ctr ) )then
                npt = npt + 1
c find intersections of contour with cell boundaries
                delta = ( ctr - v( icorn ) ) /
     +                  ( v( icorn + 1 ) - v( icorn ) )
                if( mod( icorn, 2 ) .eq. 1 )then
                  x( npt ) = i + id( icorn + 1 )
                  y( npt ) = real( j + id( icorn ) ) +
     +               delta * real( id( icorn + 1 ) - id( icorn ) )
                else
                  x( npt ) = real( i + id( icorn + 1 ) ) +
     +               delta * real( id( icorn + 2 ) - id( icorn + 1 ) )
                  y( npt ) = j + id( icorn )
                end if
                if( npt .eq. 2 )then
c draw contour segment
                  call tr( x( 1 ), y( 1 ), xx, yy )
                  call jsmove( xx, yy )
                  call tr( x( 2 ), y( 2 ), xx, yy )
                  call jsline( xx, yy )
                  npt = 0
                end if
              end if
            end do
          end do
        end do
      end do
      return
      end
