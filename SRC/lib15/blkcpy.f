      subroutine blkcpy( x, y, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c copies n consecutive real*4 values from x(i) to y(i)
c
c calling arguments
      integer n
      real x( n ), y( n )
c
c local variable
      integer i
c
      do i = 1, n
        y( i ) = x( i )
      end do
      return
      end

      subroutine blkcpy2( x, y, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c copies n consecutive real*8 values from x(i) to y(i)
c
c calling arguments
      integer n
      real*8 x( n ), y( n )
c
c local variable
      integer i
c
      do i = 1, n
        y( i ) = x( i )
      end do
      return
      end
