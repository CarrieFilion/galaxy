      subroutine trcont( x1, y1, x2, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c auxiliary routine used only by CONTOUR to transform the coordinates of
c   a point in the countour plane
c
c calling arguments
      real x1, y1, x2, y2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real grofu
c
      call trans( x1, y1, x2, y2 )
c convert to external units
      if( danl .or. dnst .or. s3dc )x2 = grofu( x2 ) / lscale
      return
      end
