      subroutine matrix( nmodes, omega, c, mdc )
c  Copyright (C) 2015, Jerry Sellwood
c
c part of mode fitting software - called from LSFUN2 and from ALPBET
c
c This routine determines the best choices of eigenvectors for the given
c   set of frequencies.  The set of linear equations with multiple RHSs
c   that determines the best choices for the eigenfunctions is derived in
c   /home/sellwood/docs/modefit.tex
c A singular matrix, or other failure, is trapped and a flag is set to
c   abort the fitting process
c   uses LAPACK routine DGESV
      use aarrays
      implicit none
c
c calling arguments
      integer nmodes, mdc
      real*8 omega( 2 * nmodes ), c( 2, mdc, * )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
c local arrays
      integer mts, mps
      parameter ( mts = 2001, mps = 201 )
      integer ipiv( 2 * mmodes )
      real*8 a( 2 * mmodes, 2 * mmodes ), b( 2 * mmodes, mps )
c
      real*8 wc, ws, wa
      common / sincos / wc( mts, mmodes ), ws( mts, mmodes ), wa( mts )
c
c local variables
      integer ifail, im, ip, it, jj, jm, kk, lp, mp, nm2, ntimes
      real*8 amp, c1, c2, phase, t
c
c check space
      nm2 = nmodes * 2
      ntimes = kt - jt + 1
      mp = kp - jp + 1
      if( nmodes .gt. mmodes )
     +                    call space( mmodes, nmodes, 'work', 'MATRIX' )
      if( ntimes .gt. mts )call space( mts, ntimes, 'work', 'MATRIX' )
      if( mp .gt. mps )call space( mp, mps, 'work', 'MATRIX' )
c compute and store cosine, sine and exponential factors
      kk = 0
      do im = 1, nmodes
        do it = 1, ntimes
          t = tme( jt + it - 1 ) - tme( kt )
          if( im .eq. 1 )wa( it ) = exp( -apodc * t )
          if( gronly )then
            phase = 0
            amp = exp( omega( kk + 1 ) * t )
          else
            phase = omega( kk + 1 ) * t
            if( ptonly )then
              amp = 1
            else
              amp = exp( omega( kk + 2 ) * t )
            end if
          end if
          wc( it, im ) = amp * cos( phase )
          ws( it, im ) = amp * sin( phase )
        end do
        kk = kk + 1
        if( .not. ( ptonly .or. gronly ) )kk = kk + 1
      end do
c set matrix of coeffs
      kk = -2
c work over rows
      do im = 1, nmodes
        kk = kk + 2
c work over columns
        jj = -2
        do jm = 1, nmodes
          jj = jj + 2
c sum over times
          c1 = 0.
          c2 = 0.
          do it = 1, ntimes
            c1 = c1 + wc( it, im ) * wc( it, jm )
     +              + ws( it, im ) * ws( it, jm )
            c2 = c2 + ws( it, im ) * wc( it, jm )
     +              - wc( it, im ) * ws( it, jm )
          end do
c set coefficients
          a( kk + 1, jj + 1 ) = c1
          a( kk + 1, jj + 2 ) = c2
          a( kk + 2, jj + 1 ) = c2
          a( kk + 2, jj + 2 ) = -c1
        end do
      end do
c set vectors of data
      do ip = 1, mp
        do im = 1, nmodes
          kk = 2 * ( im - 1 )
          b( kk + 1, ip ) = 0.
          b( kk + 2, ip ) = 0.
          lp = ip - mp
          do it = 1, ntimes
            lp = lp + mp
            b( kk + 1, ip ) = b( kk + 1, ip )
     +                    + wa( it ) * ( sdata( 1, lp ) * wc( it, im )
     +                                 + sdata( 2, lp ) * ws( it, im ) )
            b( kk + 2, ip ) = b( kk + 2, ip )
     +                    + wa( it ) * ( sdata( 1, lp ) * ws( it, im )
     +                                 - sdata( 2, lp ) * wc( it, im ) )
          end do
        end do
      end do
c solve matrix - LAPACK routine
      call dgesv( nm2, mp, a, 2 * mmodes, ipiv, b, 2 * mmodes, ifail )
c check for a failure condition
      if( ifail .ne. 0 )then
        iabort = -1
        print *, 'LSFUN2 failed for frequency vector'
        print '( 10f10.3 )', omega
        print *, 'Matrix singular for coefficients:'
        do jm = 1, nm2
          print '( 10e12.4 )', ( a( jm, im ), im = 1, nm2 )
        end do
      end if
      do im = 1, nmodes
        kk = 2 * ( im - 1 )
        do ip = 1, mp
          c( 1, im, ip ) = b( kk + 1, ip )
          c( 2, im, ip ) = b( kk + 2, ip )
        end do
      end do
      return
      end
