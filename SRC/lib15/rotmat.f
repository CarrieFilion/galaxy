      subroutine rotmat( x, e )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to rotate both position and velocity vectors using a pre-calculated
c   rotation matrix
c
c calling arguments
      real x( 6 )
      real*8 e( 3, 3 )
c
c local array
      real uc( 6 )
c
c local variables
      integer i, j
c
c save old coordinates
      do i = 1, 6
        uc( i ) = x( i )
      end do
c matrix multiplication - the first index is for the column of the matrix
      do i = 1, 3
        x( i ) = 0
        x( i + 3 ) = 0
        do j = 1, 3
          x( i ) = x( i ) + e( i, j ) * uc( j )
          x( i + 3 ) = x( i + 3 ) + e( i, j ) * uc( j + 3 )
        end do
      end do
      return
      end
