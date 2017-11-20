      subroutine propagate_parallel(x, y, m, n)
!  Copyright (C) 2014, Richard James
      
!     Purpose:

!        To propagate charges held on the edges of a parallel pair of 
!        boundary planes of the mesh across the mesh planes.  The
!        propagation direction is parallel to the store direction.        

      use converters, only: sines
      
      use interfac, local => propagate_parallel

!     Parameters:

!     x  -  a double precision rank 1 array containing values for the 
!           boundary-adjacent planes to which the charges should be
!           added.  The data is arranged as m rows of n contiguous
!           elements in memory.  Values are in transform space and all
!           interior values in the array x need modification.

!     y  -  a double precision rank 2 array containing the edge charges.
!           The second index specifies the first of the lines of edge
!           charges for propagation.

!     m  -  (integer) the number of rows in x.

!     n  -  (integer) the number of columns in x.

!     Parameter definitions:

      integer, intent(in)             :: m, n
      double precision, intent(inout) :: x(m*n)
      double precision, intent(in)    :: y(m,2)
      
!     Local variables:

      integer                         :: i, j, q, r
      double precision, allocatable   :: edge(:)

!     Acquire working memory.

      allocate(edge(m))
      r = (n - 1)/2
      q = n*m
      
!     Differences to even values of the Fourier series index.
      
      edge	  = y(:,1) - y(:,2)
      do j = 2, r
         x(j:q:n) = edge*sines(j) + x(j:q:n)
      end do
      
!     Sums to odd values of the Fourier series index.

      edge	  = y(:,1) + y(:,2)
      do j = r+1, n - 1
         x(j:q:n) = edge*sines(j) + x(j:q:n)
      end do
      
!     Release working memory.

      deallocate(edge)

!     Finish.
      
      return
      end
