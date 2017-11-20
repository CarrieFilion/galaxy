      real*8 function sphbsj( x, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c spherical Bessel functions of arbitrary order
c
c calling arguments
      integer n
      real*8 x
c
c local variables
      integer i
      real*8 ja, jb
c
c special value
      if( x .eq. 0.d0 )then
        sphbsj = 0
        if( n .eq. 0 )sphbsj = 1
      else
c first two functions
        if( n .eq. 0 )then
          sphbsj = sin( x ) / x
        else
          jb = ( sin( x ) - x * cos( x ) ) / x**2
          if( n .eq. 1 )then
            sphbsj = jb
          else
c small arguments
            if( x**n .lt. 0.1d0 )then
              i = 0
              sphbsj = 1
              do while ( i .lt. n )
                i = i + 1
                sphbsj = sphbsj * x / dble( 2 * i + 1 )
              end do
            else
c recurrence relation - it is unstable for n > x (see Press et al p 142)!
              i = 1
              ja = sin( x ) / x
              do while ( i .lt. n )
                sphbsj = dble( 2 * i + 1 ) * jb / x - ja
                i = i + 1
                ja = jb
                jb = sphbsj
              end do
            end if
          end if
        end if
      end if
      return
      end
