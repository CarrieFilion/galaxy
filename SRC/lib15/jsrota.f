      subroutine jsrota( deg )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to rotate the orientation of the text to be output subsequently
c   the angle is in counter-clockwise degrees from the x-axis
c
c calling argument
      real deg
c
c common block
c
      include 'inc/jscmmn.f'
c
      angle = deg
      return
      end
