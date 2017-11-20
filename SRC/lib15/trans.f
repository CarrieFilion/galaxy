      subroutine trans( x1, y1, x2, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine needed by contouring software to convert input position (x1,y1)
c   in units of the array indices to world coordinates (x2,y2)
c
c calling arguments
      real x1, x2, y1, y2
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      x2 = x1 * xscale + xoffs
      y2 = y1 * yscale + yoffs
      return
      end
