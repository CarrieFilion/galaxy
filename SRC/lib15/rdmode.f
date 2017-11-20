      subroutine rdmode( imode, m, renorm )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c plots the radial profile of the mode - part of mode fitting software
c
c The routine first reconstructs the mode shape from the appropriate weighted
c   sum of the basis functions and then rotates the mode so that its phase is
c   as close to zero as possible where the amplitude is high.
c If the calling argument renorm is .true., the amplitudes are rescaled so
c   that the peak is unity.
c If the function has significant phase variations, the routine draws both
c   the real (full drawn) and imaginary (dashed) parts of the function as
c   well as the enveloping amplitude (dotted).  Otherwise, only the real
c   part is drawn.
c
c calling argument
      integer imode, m
      logical renorm
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      common / sphlmn / l, mm, nn, mode
      integer l, mm, nn, mode
c
c externals
      real grofu, roundup
c
c local arrays
      real a( 101 ), ph( 101 ), r( 101 ), re( 101 ), im( 101 )
c
c local variables
      integer ir, j, k, n
      logical cplx
      real ampf, cp, p, phm, rmx, sp, w, x, y, ymn, ymx
      include 'inc/pi.f'
c
c set scaling radial grid points
      rmx = rgrid( jgrid ) / lscale
      if( sfpc )then
c basis functions
        j = 1
        k = 99
        rmx = maxr / lscale
        do ir = j, k
          r( ir ) = .01 * real( ir ) * rmx
        end do
      else if( danl .or. zanl )then
c density or bending analysis
        j = jp
        k = kp
        do ir = j, k
          if( danl )r( ir ) = grofu( real( ir ) ) / lscale
          if( zanl )r( ir ) = ( real( ir ) - .5 ) / ( drfac * lscale )
        end do
      else if( lgsp )then
c logarithmic spirals
        j = 1
        k = 61
        do ir = j, k
          r( ir ) = .01 * rmx * real( ir )
        end do
      else if( sphb )then
c spherical Bessel functions
        j = 1
        k = 50
        do ir = j, k
          r( ir ) = .02 * real( ir ) * rbess
        end do
      else
        call crash( 'RDMODE', 'Unrecognized data type' )
      end if
c get mode shape
      call mshape( imode, m, j, k, r, re, im )
c find peak amplitude
      ymx = 0
      x = 0
      do ir = j, k
        a( ir ) = sqrt( re( ir )**2 + im( ir )**2 )
        ph( ir ) = atan2( im( ir ), re( ir ) )
        if( ir .gt. 2 * m )then
          ymx = max( ymx, abs( a( ir ) ) )
          if( r( ir ) * a( ir )**2 .gt. x )then
            x = r( ir ) * a( ir )**2
            phm = ph( ir )
          end if
        end if
      end do
c find principal axis of mode
      p = 0
      ampf = 0
      do ir = j, k
        w = r( ir ) * a( ir )**2
        ampf = ampf + w
        if( ph( ir ) - phm .gt. 1.5 * pi )ph( ir ) = ph( ir ) - 2. * pi
        if( phm - ph( ir ) .gt. 1.5 * pi )ph( ir ) = ph( ir ) + 2. * pi
        if( abs( ph( ir ) - phm ) .lt. .5 * pi )then
          p = p + w * ph( ir )
        else
          if( ph ( ir ) .gt. phm )then
            p = p + w * ( ph( ir ) - pi )
          else
            p = p + w * ( ph( ir ) + pi )
          end if
        end if
      end do
      p = p / ampf
c rotate mode
      cp = cos( p )
      sp = sin( p )
      do ir = j, k
        p = re( ir ) * cp + im( ir ) * sp
        im( ir ) = im( ir  ) * cp - re( ir ) * sp
        re( ir ) = p
      end do
c determine whether imaginary part is significant
      x = 0
      y = 0
      do ir = j, k
        x = x + re( ir )**2
        y = y + im( ir )**2
      end do
c criterion compares sums of squares so amplitude difference has to be > 1000
      cplx = y .gt. 1.e-6 * x
      if( cplx )then
        ymn = -ymx
      else
c find minimum
        ymn = 0
        do ir = j, k
          ymn = min( ymn, re( ir ) )
        end do
        if( ymn .lt. 0. )ymn = min( ymn, -.1 * ymx )
      end if
      if( .not. renorm )then
        ymx = roundup( ymx )
        ymn = -roundup( -ymn )
      else
c rescale and limit data
        do ir = j, k
          a( ir ) = a( ir ) / ymx
          re( ir ) = re( ir ) / ymx
          im( ir ) = im( ir ) / ymx
c
          a( ir ) = min( a( ir ), 1.1 )
          re( ir ) = min( re( ir ), 1.1 )
          re( ir ) = max( re( ir ), -1.1 )
          im( ir ) = min( im( ir ), 1.1 )
          im( ir ) = max( im( ir ), -1.1 )
        end do
        ymn = 1.1 * ymn / ymx
        ymx = 1.1
      end if
c start new frame and draw axes
      call jspage
      if( cplx )ymn = -ymx
      call jscale( 0., rmx, ymn, ymx )
      call jsaxis( 'x', 'Radius', 1 )
      if( danl .or. lgsp .or. sfpc )call jsaxis( 'y', 'Overdensity', 1 )
      if( sphb )call jsaxis( 'y', 'Density', 1 )
      if( zanl )call jsaxis( 'y', 'z-displacement', 1 )
c write header
      call jsbldt( 'm =' )
      call jsbldi( m, 1 )
      if( lgsp )call jsbldt( 'in-plane' )
      if( zanl )call jsbldt( 'bending' )
      call jsbldt( 'mode from run no' )
      call jsbldi( irun, 5 )
      call jswrit( .1 * rmx, 1.05 * ymx )
c output eigenfrequency
      call jsbldt( 'Eigenfreq =' )
      if( .not. gronly )then
        if( abs( freq( 1, imode ) ) .gt. 5.e-4 )then
          call jsbldf( freq( 1, imode ), 5, 3 )
        else
          call jsbldi( 0, 1 )
        end if
        if( .not. ptonly )then
          if( freq( 2, imode ) .lt. 0. )then
            call jsbldt( '-' )
          else
            call jsbldt( '+' )
          end if
        end if
      end if
      if( .not. ptonly )then
        if( abs( freq( 2, imode ) ) .gt. 5.e-4 )then
          if( gronly )then
            call jsbldf( freq( 2, imode ), 5, 3 )
          else
            call jsbldf( abs( freq( 2, imode ) ), 5, 3 )
          end if
        else
          call jsbldi( 0, 1 )
        end if
        call jsbldt( 'i' )
      end if
      if( .not. cplx )call jsbldt( '- function purely real' )
      call jswrit( .2 * rmx, .06 * ymn + .94 * ymx )
c draw data
      n = k - j + 1
      call jsjoin( r( j ), re( j ), n )
      if( cplx )then
        call jsdash( 2, 2, 2, 2 )
        call jsjoin( r( j ), im( j ), n )
        call jsdash( 0, 1, 0, 1 )
        call jsjoin( r( j ), a( j ), n )
        do ir = j, k
          a( ir ) = -a( ir )
        end do
        call jsjoin( r( j ), a( j ), n )
        call jsdash( 0, 0, 0, 0 )
      end if
      return
      end
