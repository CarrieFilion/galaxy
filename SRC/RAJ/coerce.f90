      subroutine coerce(x, shift)
!  Copyright (C) 2014, Richard James
!  with changes made by Jerry Sellwood 2015
!     deleted: use timers
!     deleted calls to timer
      
!     Purpose:

!        To perform a circular shift downwards in a one dimensional array
!        in cases where the cshift function fails owing to the array size
!        exceeding its limit.

!     Parameters:

!     x      -  a double precision rank 1 array containing the data for
!               shifting.

!     shift  -  the distance by which to shift the data.

      use interfac, local => coerce

!     Parameter definitions:

      integer, intent(in)              :: shift
      integer, save                    :: number = 0
      double precision, intent(inout)  :: x(:)
      
!     Local variables:

      integer                          :: i, j, siz
      double precision, allocatable    :: duplicate(:)
      integer                          :: id
      
!     Acquire working memory.

      allocate(duplicate(shift))
      
!     Preserve the data at the beginning of the area.

      duplicate = x(1:shift)
      
!     Find the size of the array.

      siz = size(x)
      
!     Perform the downward shift.

      do i = 1, siz-shift, shift
         x(i:i+shift-1) = x(i+shift:i+2*shift-1)
      end do
      
!     Restore the preserved data at the end of the array.

      x(siz-shift+1:siz) = duplicate
      
!     Release the working memory.

      deallocate(duplicate)
      
!     Finish.

      return
      end
