      subroutine jswrit( xa, ya )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to output the current characters in the output buffer to the
c    plotting device at the point (x,y) specified in world coordinates
c
c calling arguments
      real xa, ya
c
c common block
c
      include 'inc/jscmmn.f'
c
c local variables
      character buf*10, greek*1, out*200
      integer i, j, k, n
      real h, rad, s, xs, ys
      include 'inc/pi.f'
c
c check length of string
      nchar = nchar - 1
      if( nchar .le. 0 )call crash( 'JSWRIT', 'No string to output' )
c allow enough room for text
      rad = angle * pi / 180.
      s = .7 * height * real( nchar )
      h = 2 * height
      xs = max( xa, ux1 + h * sin( rad ) / xfac )
      xs = min( xs, ux2 - ( s * cos( rad ) + h * sin( rad ) ) / xfac )
      ys = max( ya, uy1 + h * cos( rad ) / yfac )
      ys = min( ys, uy2 - ( s * sin( rad ) - h * cos( rad ) ) / yfac )
c
c copy and parse string
      i = 1
      n = 0
      do while ( i .le. nchar )
        if( ( text( i ) .eq. '_' ) .or. ( text( i ) .eq. '^' ) )then
c shift to sub- or super-script
          out( n+1:n+1 ) = char( 92 )
          out( n+2:n+2 ) = text( i )
          n = n + 2
c bracketed string
          if( text( i + 1 ) .eq. '{' )then
            k = 2
            j = i + 2
            do while ( text( j ) .ne. '}' )
              k = k + 1
              n = n + 1
              out( n:n ) = text( j )
              j = j + 1
              if( j .gt. nchar )call
     +              crash( 'JSWRIT', 'sub-string terminator not found' )
            end do
c Greek character - count the characters only
          else
            if( text( i + 1 ) .eq. char( 92 ) )then
              k = min( 10, nchar - i )
              do j = 1, k
                buf( j:j ) = text( i + 1 + j )
              end do
              call jstrns( buf, greek, k )
              k = k + 1
            else
c assume single character
              k = 1
            end if
c copy parsed string
            do j = i + 1, i + k
              n = n + 1
              out( n:n ) = text( j )
            end do
          end if
c undo shift
          out( n+1:n+1 ) = char( 92 )
          if( text( i ) .eq. '_' )then
            out( n+2:n+2 ) = '^'
          else
            out( n+2:n+2 ) = '_'
          end if
          i = i + k + 1
          n = n + 2
        else
          n = n + 1
          if( n .gt. 200 )call crash( 'JSWRIT', 'String too long' )
          out( n:n ) = text( i )
          i = i + 1
        end if
      end do
c output string
      call jsouts( out, n, xs, ys )
      nchar = 1
      return
      end
