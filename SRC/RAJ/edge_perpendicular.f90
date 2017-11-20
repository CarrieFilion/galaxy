      subroutine edge_perpendicular(x, y, m, n)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
      
!     Purpose:

!        To calculate the edge-adjacent terms of a pair of two-dimensional
!        sine transforms held in the array x.  The edges calculated are
!        those perpendicular to the store direction.

!     Parameters:

!     x       -  a two dimensional array with two rows with the transforms
!                embedded in them.  The data for each transform is held as
!                m sets of n consecutive elements.

!     y       -  a two dimensional array with four rows to receive the
!                results.  Data for the left edges is return in the rows
!                1, 3 of this array and data for the right edges in the
!                rows 2 and 4.

!     m       -  the number of rows of data in the transforms, assumed odd.

!     n       -  the number of columns of data in the transforms, assumed odd.

      use constants
      
      use converters
      
      use interfac, local => edge_perpendicular

!     Parameter definitions:
      
      integer, intent(in)           :: m, n
      double precision, intent(in)  :: x(m*n,2)
      double precision, intent(out) :: y(m, 4)
      
!     Local variables:
      
      integer                       :: i, j, ip, change, mn, ptr
      integer, save                 :: number = 0
      double precision              :: norm
      double precision, allocatable :: line(:)
      
!     Acquire working memory.

      allocate(line(m))

!     Note the size of the source arrays.

      mn = m*n

!     Clear the results area.
      
      y = 0.0d0
      
!     Set the change indicator.

      change = (n + 1)/2
      
!     Cycle over the two transforms.

      do i = 1, 2
      
         ip	= 0
         ptr	= 2*i

!        Cycle over the rows of the current transform, accumulating the
!        results in a pair of rows of y.
      
         do j = 2, n - 1
            if(j.eq.change) ptr = 2*i - 1
            ip  		= ip + 1
            y(:,ptr)		= y(:,ptr) + sines(j)*x(ip+1:ip+mn:n,i)	
         end do
	 
!        Evaluate twice the edge-adjacent values from the sums and differences.
      
         ptr        = 2*i - 1
	 line	    = y(:,ptr) - y(:,ptr+1)
         y(:,ptr)   = y(:,ptr) + y(:,ptr+1)
         y(:,ptr+1) = line
      
      end do
      
!     Normalise the results.  The normalising factor is 1/2*n.  One 
!     factor 2 in the denominator removes the factor 2 in the results
!     obtained above.  Another factor 2 remove the factor 2 in the original
!     sine transforms.  This restores the standard 2/n factor for the
!     normalisation.

      norm = 0.5d0/real(n - 1)
      y    = norm*y
      
!     Release working memory.

      deallocate(line)

!     Finish.
      
      return
      end
