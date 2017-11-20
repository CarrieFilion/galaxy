      subroutine jspenu( x, y )
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
c pgplot
      call pgmove( x, y )
c save new pen position
      xsl = x
      ysl = y
      return
      end
