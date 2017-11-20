      subroutine incare( jst )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to complete a full step for particles that may be moving with
c   fractional time steps
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
c local variables
      integer i, is, j, jz
      real asave( 3 ), csave( 6, 2 )
c
c work over particles one at a time
      do is = 1, jst
c preserve coords of first particle
        if( is .eq. 2 )then
          do i = 1, ncoor
            csave( i, 1 ) = oldc( i, 1 )
            csave( i, 2 ) = newc( i, 1 )
          end do
          do i = 1, ndimen
            asave( i ) = acc( i, 1 )
          end do
          jz = iz( 1 )
        end if
c pick out particles in guard zones
        if( iz( is ) .le. 0 )then
          if( is .gt. 1 )then
            do i = 1, ncoor
              newc( i, 1 ) = newc( i, is )
            end do
            iz( 1 ) = iz( is )
          end if
          j = 1 - iz( 1 )
          do isubst = 2, nsub( j )
c revise old coordinates
            do i = 1, ncoor
              oldc( i, 1 ) = newc( i, 1 )
            end do
c get accelerations
            call getacc( 1 )
c move forward one sub-step
            call stpgrp( 1 )
          end do
c update this particle
          if( is .gt. 1 )then
            do i = 1, ncoor
              oldc( i, is ) = oldc( i, 1 )
              newc( i, is ) = newc( i, 1 )
            end do
            do i = 1, ndimen
              acc( i, is ) = acc( i, 1 )
            end do
          end if
c end if in guard zone
        end if
c end main loop
      end do
c restore particle 1
      if( jst .gt. 1 )then
        do i = 1, ncoor
          oldc( i, 1 ) = csave( i, 1 )
          newc( i, 1 ) = csave( i, 2 )
        end do
        do i = 1, ndimen
          acc( i, 1 ) = asave( i )
        end do
        iz( 1 ) = jz
      end if
      return
      end
