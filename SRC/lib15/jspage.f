      subroutine jspage
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to advance the plotting area to the next part of the
c    device window, or to a new page if there are no unused parts
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local variables
      character flg*2
      real xr1, xr2, xt1, xt2, yr1, yr2, yt1, yt2
c
      ipic = ipic + 1
      if( ipic .gt. 1 )then
c start new page if required
        if( ipic .gt. npic )then
          ipic = 1
          if( screen )then
            call jspnta( 0., 0., 1 )
            call jsebuf
            print *, 'Press return when ready for next page'
            read '( a )', flg
          end if
c pgplot
          call pgadvance
c restore old line weight
          call jsthik( lwght )
        end if
c mark lower left corner of paper
c        call jspnta( 0., 0., 1 )
      end if
c preserve window edges
      xt1 = ( tx1 - fx1 ) / xdim
      xt2 = ( tx2 - fx1 ) / xdim
      yt1 = ( ty1 - fy1 ) / ydim
      yt2 = ( ty2 - fy1 ) / ydim
      xr1 = ( rx1 - fx1 ) / xdim
      xr2 = ( rx2 - fx1 ) / xdim
      yr1 = ( ry1 - fy1 ) / ydim
      yr2 = ( ry2 - fy1 ) / ydim
c set coordinates of window
      fx1 = mod( ipic - 1, nx )
      fx1 = vx1 + fx1 * xdim
      fx2 = fx1 + xdim
      fy1 = ny - 1 - mod( ipic - 1, npic ) / nx
      fy1 = vy1 + fy1 * ydim
      fy2 = fy1 + ydim
c update window edges
      tx1 = xt1 * fx2 + ( 1. - xt1 ) * fx1
      tx2 = xt2 * fx2 + ( 1. - xt2 ) * fx1
      ty1 = yt1 * fy2 + ( 1. - yt1 ) * fy1
      ty2 = yt2 * fy2 + ( 1. - yt2 ) * fy1
      rx1 = xr1 * fx2 + ( 1. - xr1 ) * fx1
      rx2 = xr2 * fx2 + ( 1. - xr2 ) * fx1
      ry1 = yr1 * fy2 + ( 1. - yr1 ) * fy1
      ry2 = yr2 * fy2 + ( 1. - yr2 ) * fy1
c pgplot dimensions are inches
      call pgvsiz( rx1 / 2.54, rx2 / 2.54, ry1 / 2.54, ry2 / 2.54 )
c mark upper right corner of page used so far
c      y = real( ny ) * ydim
c      call jspnta( fx2, y, 1 )
c write pgplot identifier
      if( ipic .eq. 1 )call pgiden
c re-start dashed line if required
      if( dash )then
        partd = 0
        iseg = 1
      end if
      return
c$$$c, str*2
c$$$c, y
c$$$c
c$$$      data str / ' [?38;h' /
c$$$      str(1:1) = char( 27 )
c$$$      str(2:2) = char( 12 )
c$$$c clear cifer t5 terminal screen
c$$$            print *, str
c$$$c sgs package
c$$$          call sgs_clrz
c$$$c plot10 package
c$$$      call newpag
c$$$c calcomp package
c$$$      call begplt( 99, 1 )
      end
