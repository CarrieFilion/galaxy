      subroutine lsfun2( nm2, xc, nf, fc, uiparm, urparm, ufparm )
c  Copyright (C) 2015, Jerry Sellwood
c
c part of mode fitting software
c
c Main routine called by SUMSL to determine the sum of the residuals (fc)
c   for a given input vector (xc) of eigenfrequencies.  It also calculates
c   the vector of first derivatives (gc) at this point.
c It sets up and calls subroutine MATRIX to determine the best choices of
c   eigenvectors for the given set of frequencies, and then sums the
c   residuals and derivatives by calls to RESMD2
c There are options to abort the fitting process if any of the growth rates
c   become excessive (unlikely to be physical and causes numerical problems)
c   or if the solution of the linear system fails (usually when the matrix
c   is singular).
      use aarrays
      implicit none
c
c calling arguments
      external ufparm
      integer nf, nm2, uiparm
      real*8 xc( nm2 ), fc, urparm
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
c external
      real*8 resmd2
c
c local array
      integer mps
      parameter ( mps = 201 )
      real*8 c( 2, mmodes, mps ), gc( 2 * mmodes )
c
c local variables
      integer ip, mp, nmodes
      real efac
c
c check flag
      if( iabort .lt. 0 )then
        print *,' iabort =', iabort ,' at start of LSFUN2'
        uiparm = -1
        return
      end if
c number of modes
      if( ptonly )then
        nmodes = nm2
      else
        if( gronly )then
          nmodes = nm2
        else
          nmodes = nm2 / 2
        end if
c check for excessive growth rates
        fc = 1.d20
        do mp = 1, nmodes
          ip = mp
          if( .not. gronly )ip = 2 * mp
          efac = xc( ip ) * ( tme( kt ) - tme( jt ) )
          if( efac .gt. 20. )then
c transmit flag to abort fitting process
c$$$            print *, 'Excessive growth factor for mode no', mp
c$$$            print *, 'growth rate =', xc( ip ), ' exponent =', efac
c$$$            print *, 'Fit stopped in LSFUN2'
            uiparm = -1
            return
          end if
        end do
      end if
c check space
      if( nmodes .gt. mmodes )call space(
     +                                mmodes, nmodes, 'work', 'lsfun2' )
      mp = kp - jp + 1
      if( mp .gt. mps )call space( mps, mp, 'work', 'lsfun2' )
c find best choices of initial amplitudes and phases
      call matrix( nmodes, xc, c, mmodes )
c determine sum of residuals and gradients
      do ip = 1, nm2
        gc( ip ) = 0.
      end do
      fc = 0.
      do ip = 1, mp
        fc = fc + resmd2( uiparm, nmodes, c( 1, 1, ip ), ip, gc )
      end do
      return
      end

      subroutine lsgrd2( nm2, xc, nf, gc, uiparm, urparm, ufparm )
c  Copyright (C) 2015, Jerry Sellwood
c
c part of mode fitting software
c
      use aarrays
      implicit none
c
c calling arguments
      external ufparm
      integer nf, nm2, uiparm
      real*8 xc( nm2 ), gc( nm2 ), urparm
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
c external
      real*8 resmd2
c
c local array
      integer mps
      parameter ( mps = 201 )
      real*8 c( 2, mmodes, mps )
c
c local variables
      integer ip, mp, nmodes
      real efac
      real*8 fc
c
c check flag
      if( iabort .lt. 0 )then
        print *,' iabort =', iabort ,' at start of LSGRD2'
        uiparm = -1
        return
      end if
c number of modes
      if( ptonly )then
        nmodes = nm2
      else
        if( gronly )then
          nmodes = nm2
        else
          nmodes = nm2 / 2
        end if
c check for excessive growth rates
        do mp = 1, nmodes
          ip = mp
          if( .not. gronly )ip = 2 * mp
          efac = xc( ip ) * ( tme( kt ) - tme( jt ) )
          if( efac .gt. 20. )then
c transmit flag to abort fitting process
            print *, 'Excessive growth factor for mode no', mp
            print *, 'growth rate =', xc( ip ), ' exponent =', efac
            print *, 'Fit stopped in LSFUN2'
            uiparm = -1
            return
          end if
        end do
      end if
c check space
      if( nmodes .gt. mmodes )call space(
     +                                mmodes, nmodes, 'work', 'lsgrd2' )
      mp = kp - jp + 1
      if( mp .gt. mps )call space( mps, mp, 'work', 'lsfun2' )
c find best choices of initial amplitudes and phases
      call matrix( nmodes, xc, c, mmodes )
c transmit error flag
      if( in .lt. 0 )then
        uiparm = -1
        return
      end if
c determine sum of residuals and gradients
      do ip = 1, nm2
        gc( ip ) = 0.
      end do
      fc = 0.
      do ip = 1, mp
        fc = fc + resmd2( uiparm, nmodes, c( 1, 1, ip ), ip, gc )
      end do
      return
      end
