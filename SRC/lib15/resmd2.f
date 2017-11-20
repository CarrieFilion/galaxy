      real*8 function resmd2( iflag, nmodes, c, ip, gc )
c  Copyright (C) 2015, Jerry Sellwood
c
c part of mode fitting software - called from LSFUN2
c
c Routine to sum the residuals, via the function value, and to determine
c   the vector of derivatives (gc) for the current set of frequencies.  The
c   best choice of eigenvectors was determined by a previous call to MATRIX
      use aarrays
      implicit none
c
c calling arguments
      integer iflag, nmodes, ip
      real*8 c( 2, nmodes ), gc( 2 * nmodes )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      integer mts
      parameter ( mts = 2001 )
      common / sincos / wc( mts, mmodes ), ws( mts, mmodes ), wa( mts )
      real*8 wc, ws, wa
c
c local variables
      integer im, it, k, lp, mp, ntimes
      real*8 eta, t, xi, xx, yy
c
c compute sum of squares
      ntimes = kt - jt + 1
      resmd2 = 0.
      mp = kp - jp + 1
      lp = ip - mp
c find apodised data
      do it = 1, ntimes
        lp = lp + mp
        xx = wa( it ) * sdata( 1, lp )
        yy = wa( it ) * sdata( 2, lp )
c subtract off modes
        do im = 1, nmodes
          xx = xx - c( 1, im ) * wc( it, im ) + c( 2, im ) * ws( it, im)
          yy = yy - c( 1, im ) * ws( it, im ) - c( 2, im ) * wc( it, im)
        end do
        resmd2 = resmd2 + xx * xx + yy * yy
c compute derivatives
        t = tme( it + jt - 1 ) - tme( jt )
        k = 0
        do im = 1, nmodes
          eta = c( 1, im ) * wc( it, im ) - c( 2, im ) * ws( it, im )
          xi =  c( 1, im ) * ws( it, im ) + c( 2, im ) * wc( it, im )
          k = k + 1
          if( gronly )then
            gc( k ) = gc( k ) - 2. * t * ( xx * eta + yy * xi )
          else
            gc( k ) = gc( k ) + 2. * t * ( xx * xi - yy * eta )
            if( .not. ptonly )then
              k = k + 1
              gc( k ) = gc( k ) - 2. * t * ( xx * eta + yy * xi )
            end if
          end if
        end do
      end do
      return
      end
