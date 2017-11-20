      subroutine bcsmth( x, n, ks )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c smooths a 1-D data array using a running average
c   output data overwrite input
c
c calling arguments
      integer n, ks
      real x( n )
c
c local array
      real, allocatable :: y( : )
c
c local variables
      integer i, j, k, l1, l2
      real d
c
      if( ks .le. 0 )call crash( 'BCSMTH', 'nonsense argument' )
      k = min( ks, n )
c allocate local space
      allocate ( y( ks ) )
      do i = 1, n
c force one-sided averages at ends of range
        l1 = max( 1, i - k / 2 )
        l2 = min( i + k / 2, n )
c compute average
        d = 0
        do j = l1, l2
          d = d + x( j )
        end do
        d = d / real( l2 - l1 + 1 )
        if( i .le. ks )then
c load circular buffer
          y( i ) = d
        else
c circular buffer full
          j = mod( i - 1, ks ) + 1
          x( i - ks ) = y( j )
          y( j ) = d
        end if
      end do
c collect last results from circular buffer
      k = max( 1, n - ks + 1 )
      do i = k, n
        j = mod( i - 1, ks ) + 1
        x( i ) = y( j )
      end do
      deallocate ( y )
      return
      end
