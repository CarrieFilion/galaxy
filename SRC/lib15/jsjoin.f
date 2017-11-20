      subroutine jsjoin( x, y, n )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c primitive routine to join the points in the input arrays of x and y
c   with straight lines
c
c calling arguments
      integer n
      real x( n ), y( n )
c
c local variables
      integer i
      real xx, yy
c
      do i = 1, n
        xx = x( i )
        yy = y( i )
        if( i .eq. 1 )call jsmove( xx, yy )
        call jsline( xx, yy )
      end do
      return
      end

      subroutine jscoin( x, y, n, k )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to join the points in the input arrays of x and y with straight
c   lines that have the color specified by the input value k
c
c calling arguments
      integer n, k
      real x( n ), y( n )
c
c common block
c
      include 'inc/comprns.f'
c
      if( nruns .gt. 1 )call pgsci( k + 1 )
      call jsjoin( x, y, n )
      if( nruns .gt. 1 )call pgsci( 1 )
      return
      end
