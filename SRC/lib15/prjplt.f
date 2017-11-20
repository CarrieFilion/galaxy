      subroutine prjplt
c  Copyright (C) 2015, Jerry Sellwood
c
c Bare bones routine to contour projected density
c
c Called from ANALYS
c Graphics routines from JSPLOT
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/anlys.f'
c
c external
      external trans
c
c local variables
      integer ifail
      real xmax
c
      if( .not. c3d )call crash( 'PRJPLT', 'Not C3D grid' )
      call nextrec( 'PROJ', ifail )
      do while ( ifail .eq. 0 )
        call jspage
c x-y projection only for now
        xmax = rgrid( jgrid )
        call jsescl( -xmax, xmax, -xmax, xmax )
        call jsaxis( 'x', 'x', 1 )
        call jsaxis( 'y', 'y', 1 )
        call jsctra( 5, wres, ngx, ngy, trans, .true., .false. )
        call nextrec( 'PROJ', ifail )
      end do
      return
      end
