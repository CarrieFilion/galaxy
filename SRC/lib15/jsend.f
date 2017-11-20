      subroutine jsend
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to close the current plotting device
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local variable
      character a*1
c
      if( screen )then
        call jsebuf
        print *, 'Enter return to close plot'
        read '( a )', a
c      else
      end if
c pgplot
      call pgend
      return
      end
