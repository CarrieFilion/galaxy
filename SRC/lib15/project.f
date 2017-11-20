      subroutine project( imode, mact )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to contour the projected density of a fitted 3D mode
c   from an expansion in basis functions
c part of the mode fitting software
c
c calling argument
      integer imode, mact
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c      include 'inc/anlsis.f'
c
      include 'inc/anlys.f'
c
      common / sphlmn / l, m, n, mode
      integer l, m, n, mode
c
c externals
      external projx, projy, projz, trans
c
c local array
      integer dim
      parameter ( dim = 51 )
      real w( dim, dim )
c
c local variables
      real rmx, x1, x2, x3, x4, y1, y2, y3, y4
c
      rmx = rbess / lscale
      m = mact
      mode = imode
c
      x1 = .15
      x2 = x1 + .4
      x3 = x2
      x4 = x3 + .4
      y1 = .1
      y2 = y1 + .4
      y3 = y2
      y4 = y3 + .4
      call jssize( x3, x4, y3, y4 )
      call jscale( -1., 1., -1., 1. )
      call jsbldt( 'mode' )
      call jsbldi( imode, 1 )
      call jswrit( -.8, .5 )
c
      call jssize( x1, x2, y1, y2 )
      call jscale( -rmx, rmx, -rmx, rmx )
      call jsaxis( 'x', 'x', 1 )
      call jsaxis( 'y', 'y', 1 )
      call jscirc( 0., 0., rmx )
      call jscale( -1., 1., -1., 1. )
      call jsctrf( projz, 10, w, dim, dim, trans, .true. )
c
      call jssize( x3, x4, y1, y2 )
      call jscale( -rmx, rmx, -rmx, rmx )
      call jsaxis( 'x', 'z', 1 )
      call jsaxis( 'y', ' ', 0 )
      call jscirc( 0., 0., rmx )
      call jscale( -1., 1., -1., 1. )
      call jsctrf( projx, 10, w, dim, dim, trans, .true. )
c
      call jssize( x1, x2, y3, y4 )
      call jscale( -rmx, rmx, -rmx, rmx )
      call jsaxis( 'x', ' ', 0 )
      call jsaxis( 'y', 'z', 1 )
      call jscirc( 0., 0., rmx )
      call jscale( -1., 1., -1., 1. )
      call jsctrf( projy, 10, w, dim, dim, trans, .true. )
c
      return
      end
