      integer nx, ny, npic, ipic
      real xdim, ydim, fx1, fx2, fy1, fy2, xoffjs, yoffjs
c f for the full size of the whole, or divided part, the plotting surface
      common / jsarng / nx, ny, npic, xdim, ydim, ipic, fx1, fx2, fy1,
     +                  fy2, xoffjs, yoffjs
c
      logical printx, laserp, laserl, cifert5, tektrx, versa, pericom,
     +        graphon, args, zeta, null
      common / jsdev / printx, laserp, laserl, cifert5, tektrx, versa,
     +                 pericom, graphon, args, zeta, null
c
      integer ntab
      parameter ( ntab = 21 )
      real xtab( ntab ), ytab( ntab )
      common / jsintr / xtab, ytab
c
      real height, angle
      common / jslett / height, angle
c
      real xlast, xprev, ylast, yprev
      logical outbox
      common / jslpos / xlast, ylast, xprev, yprev, outbox
c
      integer iseg, lwght
      logical dash
      real segl( 4 ), total, partd
      common / jsltyp / lwght, dash, segl, total, partd, iseg
c
      logical screen
      common / jsmode / screen
c t for the active size of the current part of the plotting surface in cm
c r for the actual size, which could be less is equal scales are requested
c s for the bounds in world coordinates
      real ejsx, ejsy, rx1, rx2, ry1, ry2, sx1, sx2, sy1, sy2, tx1, tx2
      real ty1, ty2, xfac, yfac
      common / jspict / rx1, rx2, ry1, ry2, sx1, sx2, sy1, sy2, xfac,
     +                  yfac, tx1, tx2, ty1, ty2, ejsx, ejsy
c v for the max size of the plotting surface in cm, and u for the same in
c   world coordinates
      real ux1, ux2, uy1, uy2, vx1, vx2, vy1, vy2
      common / pgpict / vx1, vx2, vy1, vy2, ux1, ux2, uy1, uy2
c
      real pixpcm, xcmmax, ycmmax
      common / jspsiz / pixpcm, xcmmax, ycmmax
c
      real xsl, ysl
      common / jsppos / xsl, ysl
c
      integer ibuff, nchar
      parameter ( ibuff = 200 )
      character*1 text( ibuff )
      common / jstrng / text
      common / jstrng2 / nchar
c
c      common / sgspar / izonid, istat
c      integer izonid, istat
