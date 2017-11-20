      subroutine jspend( x, y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      real x, y
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local variables
      real dx, dy, rs
c
c compute line length in pixels
      dx = ( x - xsl ) * ( rx2 - rx1 ) / ( sx2 - sx1 )
      dy = ( y - ysl ) * ( ry2 - ry1 ) / ( sy2 - sy1 )
      rs = pixpcm * sqrt( dx**2 + dy**2 )
c dot required
      if( rs .lt. .1 )then
        call jspnta( x, y, lwght )
      else
c pgplot
        call pgdraw( x, y )
      end if
c save new pen position
      xsl = x
      ysl = y
      return
      end
