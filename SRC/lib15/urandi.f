      subroutine urandi( ia, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to fill the given array with randomly arranged indices in the
c   range 1 - n
c
c calling arguments
      integer n, ia( n )
c
c external
      real*8 ranuni
c
c local allocatable array
      real, allocatable :: w(:)
c
c local variables
      integer i
c
      allocate ( w( n ) )
c generate n random numbers
      do i = 1, n
        w( i ) = ranuni( 0. )
      end do
c rank them
      call rnkmrg( w, n, ia )
      deallocate ( w )
      return
      end
