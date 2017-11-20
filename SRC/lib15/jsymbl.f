      subroutine jsymbl( x, y, n, k )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to output symbols having the style specified by k at the 
c   n positions (x,y)
c
c calling arguments
      integer n, k
      real x( n ), y( n )
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local variables
      integer i, j
      real sizex, sizey, th, xx, x1, yy, y1
      include 'inc/pi.f'
c
c set half symbol size
      sizex = .5 * height * ( sx2 - sx1 ) / ( rx2 - rx1 )
      sizey = .5 * height * ( sy2 - sy1 ) / ( ry2 - ry1 )
c work over points
      do i = 1, n
        call jscoor( x( i ), y( i ), xx, yy )
c discard points outside current window
        if( ( xx .gt. sx1 ) .and. ( xx .lt. sx2 ) .and.
     +      ( yy .gt. sy1 ) .and. ( yy .lt. sy2 ) )then
          if( k .eq. 0 )then
c dot
            call jspnta( xx, yy, lwght )
          else if( k .eq. 1 )then
c square
            call jspenu( xx + sizex, yy + sizey )
            call jspend( xx - sizex, yy + sizey )
            call jspend( xx - sizex, yy - sizey )
            call jspend( xx + sizex, yy - sizey )
            call jspend( xx + sizex, yy + sizey )
          else if( k .eq. 2 )then
c circle
            do j = 1, 21
              th = .1 * pi * real( j - 1 )
              x1 = xx + sizex * cos( th )
              y1 = yy + sizey * sin( th )
              if( j .eq. 1 )call jspenu( x1, y1 )
              call jspend( x1, y1 )
            end do
          else if( k .eq. 3 )then
c triangle
            call jspenu( xx + sizex, yy - sizey )
            call jspend( xx, yy + sizey )
            call jspend( xx - sizex, yy - sizey )
            call jspend( xx + sizex, yy - sizey )
          else if( k .eq. 4 )then
c plus
            call jspenu( xx, yy - sizey )
            call jspend( xx, yy + sizey )
            call jspenu( xx - sizex, yy )
            call jspend( xx + sizex, yy )
          else if( k .eq. 5 )then
c cross
            call jspenu( xx - sizex, yy - sizey )
            call jspend( xx + sizex, yy + sizey )
            call jspenu( xx - sizex, yy + sizey )
            call jspend( xx + sizex, yy - sizey )
          else if( k .eq. 6 )then
c diamond
            call jspenu( xx, yy - sizey )
            call jspend( xx + sizex, yy )
            call jspend( xx, yy + sizey )
            call jspend( xx - sizex, yy )
            call jspend( xx, yy - sizey )
          else
            call crash( 'JSYMBL', 'Symbol index out of range' )
          end if
        end if
      end do
      return
      end
