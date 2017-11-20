      subroutine jsrscl( x1, x2, y1, y2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to rescale the input values (x1,y1) to output (x2,y2)
c   when using 2D kernel estimates of a density specified by a
c   number of particles.  It is needed to ensure that the kernel
c   is round even if the x- and y-axes are scale differently
c
c calling arguments
      real x1, x2, y1, y2
c
c common block
c
      include 'inc/jscmmn.f'
c
c local variables
      real a1, a2
c
      sx1 = x1
      sx2 = x2
      sy1 = y1
      sy2 = y2
      xfac = ( tx2 - tx1 ) / ( sx2 - sx1 )
      yfac = ( ty2 - ty1 ) / ( sy2 - sy1 )
      ejsx = 1.e-5 * ( sx2 - sx1 )
      ejsy = 1.e-5 * ( sy2 - sy1 )
c reverse engineer pgplot scaling - needed only for text writing
      a1 = ( rx1 - vx1 ) / ( vx2 - vx1 )
      a2 = ( rx2 - vx1 ) / ( vx2 - vx1 )
      ux1 = ( a2 * sx1 - a1 * sx2 ) / ( a2 - a1 )
      ux2 = ( ( a2 - 1. ) * sx1 - ( a1 - 1 ) * sx2 ) / ( a2 - a1 )
      a1 = ( ry1 - vy1 ) / ( vy2 - vy1 )
      a2 = ( ry2 - vy1 ) / ( vy2 - vy1 )
      uy1 = ( a2 * sy1 - a1 * sy2 ) / ( a2 - a1 )
      uy2 = ( ( a2 - 1. ) * sy1 - ( a1 - 1 ) * sy2 ) / ( a2 - a1 )
c      call pgswin( ux1, ux2, uy1, uy2 )
c
c pgplot dimensions are inches
      call pgvsiz( rx1 / 2.54, rx2 / 2.54, ry1 / 2.54, ry2 / 2.54 )
      call pgswin( sx1, sx2, sy1, sy2 )
      return
      end
