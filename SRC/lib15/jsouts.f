      subroutine jsouts( in, nchars, xs, ys )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to output the character string "in" to the plotting device
c   the specifications of super- and sub-scripts, and of Greek characters
c   are converted from TEX format to that required by PGPLOT
c
c calling arguments
      integer nchars
      character*(*) in
      real xs, ys
c
c common block
c
      include 'inc/jscmmn.f'
c
c local array
      character*200 out
c
c local variables
      character*1 greek
      integer i, k, n
c
c parse string
      k = 0
      i = 1
      do while ( i .le. nchars )
c backslash flags a special sequence
        if( ( i .lt. nchars ) .and. ( in( i:i ) .eq. char( 92 ) ) )then
          if( k + 3 .gt. 200 )call crash( 'JSOUTS', 'String too long' )
          out( k+1:k+1 ) = in( i:i )
c sub- or super-script indicator
          if( ( in( i+1:i+1 ) .eq. '_' ) .or.
     +        ( in( i+1:i+1 ) .eq. '^' ) )then
            if( in( i+1:i+1 ) .eq. '_' )out( k+2:k+2 ) = 'd'
            if( in( i+1:i+1 ) .eq. '^' )out( k+2:k+2 ) = 'u'
            k = k + 2
            i = i + 2
          else
c Greek character
            call jstrns( in( i+1:nchars ), greek, n )
            out( k+2:k+2 ) = 'g'
            out( k+3:k+3 ) = greek
            k = k + 3
            i = i + 1 + n
          end if
        else
c copy Roman characters
          k = k + 1
          if( k .gt. 200 )call crash( 'JSOUTS', 'String too long' )
          out( k:k ) = in( i:i )
          i = i + 1
        end if
      end do
c pgplot
      call pgptext( xs, ys, angle, 0., out( 1:k ) )
      return
      end
