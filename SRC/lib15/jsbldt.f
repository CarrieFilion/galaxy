      subroutine jsbldt( chars )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to add a character string to the current output buffer
c
c calling argument
      character*(*) chars
c
c common block
c
      include 'inc/jscmmn.f'
c
c local variables
      character a*1, bslash*2
      integer i, n
      parameter ( a = ' ', bslash = char( 92 ) // char( 92 ) )
c
c count number of characters in argument string
      n = len( chars )
c truncate string if total exceeds 200 chars
      n = min( n, ibuff - nchar )
      if( n .le. 0 )return
c copy new text into buffer
      i = 1
      do while ( i .le. n )
c eliminate any double backslashes - some compilers pass both, others don't
        if( ( i .lt. n ) .and.
     +      ( chars( i:i+1 ) .eq. bslash ) )then
          i = i + 1
        end if
        text( nchar ) = chars( i:i )
        nchar = nchar + 1
        i = i + 1
      end do
c add a blank character after string
      if( nchar .lt. ibuff )then
        text( nchar ) = a
        nchar = nchar + 1
      end if
      return
      end
