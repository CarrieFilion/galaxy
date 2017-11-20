      subroutine shuffl( a, nc, nr )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c utility routine to randomly shuffle columns of a matrix
c
c calling arguments
      integer nc, nr
      real a( nr, nc )
c
c local allocatable arrays
      integer, allocatable :: ind(:)
      real, allocatable :: b(:,:)
c
c local variables
      integer i, j, k
c
      allocate ( b( nr, nc ) )
      allocate ( ind( nc ) )
c generate an array of random, but unique, indices
      call urandi( ind, nc )
c copy each row from random columns in a to sequential columns in b
      do j = 1, nc
        k = ind( j )
        do i = 1, nr
          b( i, j ) = a( i, k )
        end do
      end do
      deallocate( ind )
c copy back
      call blkcpy( b( 1, 1 ), a( 1, 1 ), nc * nr )
      deallocate( b )
      return
      end
