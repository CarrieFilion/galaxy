      real*8 function ranuni( dummy )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns a random sample from a uniform distribution between 0 and 1
c   uses the NAg generator G05CAF
c
c calling argument unused
      real dummy
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real*8 g05caf
c
c local array
      integer mbuf
      parameter ( mbuf = 100 )
      real x( mbuf )
      real*8 ran( mbuf )
c
c local variables
      integer i, j
      save j, ran
      data j / mbuf /
c
      if( lnag )then
c random generator uniform in range 0 - 1
        ranuni = g05caf( dummy )
      else
c        call nonag( 'RANUNI', 'G05CAF' )
c f90 random generator
        j = j + 1
        if( j .gt. mbuf )then
c fill buffer with new values
          call random_number( x )
          do i = 1, mbuf
            ran( i ) = x( i )
          end do
          j = 1
        end if
        ranuni = ran( j )
      end if
      return
      end
