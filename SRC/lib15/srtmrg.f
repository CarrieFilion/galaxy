      subroutine srtmrg( a, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to sort a vector of real data into ascending order using a
c   local workspace.  It is highly efficient - O(NlogN)
c
c calling arguments
      integer n
      real a( n )
c
c local allocatable array
      real, allocatable :: w(:)
c
c local variables
      integer i, ip, is, j, k, l, m
c
      if( n .ge. 2 )then
        allocate ( w( n / 2 ) )
c iterative merge
        is = 1
        do while ( is .lt. n )
          i = n + 1
          do while ( i .gt. is )
c i points to the 1st element of the 1st part, j to the 1st of the 2nd part
            j = i - is
            i = max( i - 2 * is, 1 )
c always is elements in 2nd part, but there may be fewer in 1st part
            k = min( j - i, is )
            if( k .gt. 0 )then
c merge only if the values if the two parts overlap
              if( a( j - 1 ) .gt. a( j ) )then
c duplicate first part
                do l = 1, k
                  w( l ) = a( i - 1 + l )
                end do
c merge the two parts in rank order
                l = 1
                m = 1
                ip = i
c interleave first and second parts
                do while ( ( l .le. k ) .and. ( m .le. is ) )
                  if( w( l ) .le. a( j ) )then
                    a( ip ) = w( l )
                    l = l + 1
                  else
                    a( ip ) = a( j )
                    m = m + 1
                    j = j + 1
                  end if
                  ip = ip + 1
                end do
c copy any elements of w that are greater than the last element of 2nd part
                do while ( l .le. k )
                  a( ip ) = w( l )
                  l = l + 1
                  ip = ip + 1
                end do
              end if
            end if
          end do
c double vector lengths for next level
          is = 2 * is
        end do
        deallocate ( w )
      end if
      return
      end

      subroutine srtmrg2( a, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to sort a vector of real*8 data into ascending order using a
c   local workspace.  It is highly efficient - O(NlogN)
c
c calling arguments
      integer n
      real*8 a( n )
c
c local allocatable array
      real*8, allocatable :: w(:)
c
c local variables
      integer i, ip, is, j, k, l, m
c
      if( n .ge. 2 )then
        allocate ( w( n / 2 ) )
c iterative merge
        is = 1
        do while ( is .lt. n )
          i = n + 1
          do while ( i .gt. is )
c i points to the 1st element of the 1st part, j to the 1st of the 2nd part
            j = i - is
            i = max( i - 2 * is, 1 )
c always is elements in 2nd part, but there may be fewer in 1st part
            k = min( j - i, is )
            if( k .gt. 0 )then
c merge only if the values if the two parts overlap
              if( a( j - 1 ) .gt. a( j ) )then
c duplicate first part
                do l = 1, k
                  w( l ) = a( i - 1 + l )
                end do
c merge the two parts in rank order
                l = 1
                m = 1
                ip = i
c interleave first and second parts
                do while ( ( l .le. k ) .and. ( m .le. is ) )
                  if( w( l ) .le. a( j ) )then
                    a( ip ) = w( l )
                    l = l + 1
                  else
                    a( ip ) = a( j )
                    m = m + 1
                    j = j + 1
                  end if
                  ip = ip + 1
                end do
c copy any elements of w that are greater than the last element of 2nd part
                do while ( l .le. k )
                  a( ip ) = w( l )
                  l = l + 1
                  ip = ip + 1
                end do
              end if
            end if
          end do
c double vector lengths for next level
          is = 2 * is
        end do
        deallocate ( w )
      end if
      return
      end
