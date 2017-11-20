      subroutine rnkmrg( a, n, irank )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to rank a vector of real data into ascending order using a local
c   workspace.  It is highly efficient - O(NlogN)
c The input values of a are returned unchanged, while the returned value
c   of each element irank( i ) = j, where a( j ) is the j'th largest
c   value in the input array
c
c calling arguments
      integer n, irank( n )
      real a( n )
c
c local allocatable array
      integer, allocatable :: iwork(:)
c
c local variables
      integer i, ip, is, j, k, l, m
c
c initialize
      do i = 1, n
        irank( i ) = i
      end do
      if( n .ge. 2 )then
        allocate ( iwork( n / 2 ) )
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
              if( a( irank( j - 1 ) ) .gt. a( irank( j ) ) )then
c duplicate first part
                do l = 1, k
                  iwork( l ) = irank( i - 1 + l )
                end do
c merge the two parts in rank order
                l = 1
                m = 1
                ip = i
c interleave first and second parts
                do while ( ( l .le. k ) .and. ( m .le. is ) )
                  if( a( iwork( l ) ) .le. a( irank( j ) ) )then
                    irank( ip ) = iwork( l )
                    l = l + 1
                  else
                    irank( ip ) = irank( j )
                    m = m + 1
                    j = j + 1
                  end if
                  ip = ip + 1
                end do
c copy any elements of iwork that are greater than the last element of 1st part
                do while ( l .le. k )
                  irank( ip ) = iwork( l )
                  l = l + 1
                  ip = ip + 1
                end do
              end if
            end if
          end do
          is = 2 * is
        end do
        deallocate ( iwork )
      end if
      return
      end

      subroutine rnkmrg2( a, n, irank )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to rank a vector of real*8 data into ascending order using a
c   local workspace.  It is highly efficient - O(NlogN)
c The input values of a are returned unchanged, while the returned value
c   of each element irank( i ) = j, where a( j ) is the j'th largest
c   value in the input array
c
c calling arguments
      integer n, irank( n )
      real*8 a( n )
c
c local allocatable array
      integer, allocatable :: iwork(:)
c
c local variables
      integer i, ip, is, j, k, l, m
c
c initialize
      do i = 1, n
        irank( i ) = i
      end do
      if( n .ge. 2 )then
        allocate ( iwork( n / 2 ) )
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
              if( a( irank( j - 1 ) ) .gt. a( irank( j ) ) )then
c duplicate first part
                do l = 1, k
                  iwork( l ) = irank( i - 1 + l )
                end do
c merge the two parts in rank order
                l = 1
                m = 1
                ip = i
c interleave first and second parts
                do while ( ( l .le. k ) .and. ( m .le. is ) )
                  if( a( iwork( l ) ) .le. a( irank( j ) ) )then
                    irank( ip ) = iwork( l )
                    l = l + 1
                  else
                    irank( ip ) = irank( j )
                    m = m + 1
                    j = j + 1
                  end if
                  ip = ip + 1
                end do
c copy any elements of iwork that are greater than the last element of 1st part
                do while ( l .le. k )
                  irank( ip ) = iwork( l )
                  l = l + 1
                  ip = ip + 1
                end do
              end if
            end if
          end do
          is = 2 * is
        end do
        deallocate ( iwork )
      end if
      return
      end
