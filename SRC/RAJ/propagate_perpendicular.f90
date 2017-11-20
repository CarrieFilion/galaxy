      subroutine propagate_perpendicular(x, y, m, n)
!  Copyright (C) 2014, Richard James
      
!     Purpose:

!        To propagate charges held on the edges of a parallel pair of 
!        boundary planes of the mesh across the boundary planes.  The
!        popagation direction is perpendicular to the store direction.        

      use converters, only: sines
      
      use interfac, local => propagate_perpendicular
      
!     Parameters:

!     x  -  a double precision rank 1 array containing values for the 
!           boundary-adjacent planes to which the charges should be
!           added.  The data is arranged as m rows of n contiguous
!           elements in memory.  Values are in transform space and all
!           interior values of x need modification.

!     y  -  a double precision rank 2 array containing the edge charges.
!           The second index specifies the first line of charges for
!           propagation.

!     m  -  (integer) the number of rows in x.

!     n  -  (integer) the number of columns in x.

!     Parameter definitions:

      integer, intent(in)             :: m, n
      double precision, intent(inout) :: x(m*n)
      double precision, intent(in)    :: y(n,2)
      
!     Local variables:

      integer                         :: i, j, p, q, r
      integer, save                   :: number = 0
      double precision, allocatable   :: edge(:)

!     Acquire working memory.

      allocate(edge(n))
      r = (m - 1)/2
      p = n + 1 	  !	  Start of line in x.
      q = 2*n		  !	  End of line in x.
      
!     Differences to even values of the Fourier series index.
      
      edge	= y(:,1) - y(:,2)
      do j = 2, r
         x(p:q) = edge*sines(j) + x(p:q)
         p	= p + n
         q	= q + n
      end do
      
!     Sums to odd values of the Fourier series index.

      edge	= y(:,1) + y(:,2)
      do j = r+1, m - 1
         x(p:q) = edge*sines(j) + x(p:q)
         p	= p + n
         q	= q + n
      end do
      
!     Release working memory.

      deallocate(edge)

!     Finish.
      
      return
      end
