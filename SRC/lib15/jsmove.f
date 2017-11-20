      subroutine jsmove( x, y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to move the pen to the new point (x,y) without drawing a line
c
c calling arguments
      real x, y
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c move pen
      call jspenu( x, y )
c save new position
      xlast = x
      ylast = y
      if( dash )then
        partd = 0.
        iseg = 1
      end if
      return
      end
