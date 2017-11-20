      subroutine alpbet( nmodes, sumsq )
c  Copyright (C) 2015, Jerry Sellwood
c
c part of mode fitting software
c
c routine to re-determine, and save in useful form, the best fitting complex
c   eigenfunctions for the final optimal frequencies
c it also recomputes the sum of the squared residuals which needs to be
c   re-calculated if the data were re-scaled before fitting
      use aarrays
      implicit none
c
c calling arguments
      integer nmodes
      real sumsq
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
c local arrays
      integer mps
      parameter ( mps = 201 )
      real*8 omega( 2 * mmodes ), c( 2, mmodes, mps )
c
c local variables
      integer imode, ip, it, k, lp, mp, n
      real cp, efac, phase, sp, t, x, y
c
c check space
      if( nmodes .gt. mmodes )
     +                    call space( mmodes, nmodes, 'work', 'ALPBET' )
      mp = kp - jp + 1
      if( mp .gt. mps )call space( mps, mp, 'work', 'ALPBET' )
c set frequencies
      k = 0
      do imode = 1, nmodes
        k = k + 1
        if( gronly )then
          omega( k ) = freq( 2, imode ) - apodc
        else
          omega( k ) = freq( 1, imode )
          if( .not. ptonly )then
            k = k + 1
            omega( k ) = freq( 2, imode ) - apodc
          end if
        end if
      end do
c find best choices of initial amplitudes and phases
      call matrix( nmodes, omega, c, mmodes )
c save these
      do lp = 1, mp
        ip = lp + jp - 1
        do imode = 1, nmodes
          alp( ip, imode ) = c( 1, imode, lp )
          bet( ip, imode ) = c( 2, imode, lp )
        end do
      end do
c recompute sumsq if apodc is non-zero
      if( apodc .ne. 0. )then
        sumsq = 0.
        n = 0
        do it = jt, kt
          t = tme( it ) - tme( kt )
          do ip = jp, kp
c compute fitted value
            x = 0
            y = 0
            do imode = 1, nmodes
              phase = t * freq( 1, imode )
              cp = cos( phase )
              sp = sin( phase )
              efac = exp( t * freq( 2, imode ) )
              x = x + efac * ( alp( ip, imode ) * cp
     +                       - bet( ip, imode ) * sp )
              y = y + efac * ( alp( ip, imode ) * sp
     +                       + bet( ip, imode ) * cp )
            end do
c compute residual
            n = n + 1
            x = x - sdata( 1, n )
            y = y - sdata( 2, n )
            sumsq = sumsq + x * x + y * y
          end do
        end do
        sumsq = sumsq / denf
      end if
      print '( '' Normalised sum of squares ='', f8.4 )', sumsq
      return
      end
