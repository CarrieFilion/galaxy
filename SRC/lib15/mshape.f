      subroutine mshape( imode, m, j, k, r, re, im )
c  Copyright (C) 2015, Jerry Sellwood
c
c determines the radial profile of the mode - part of mode fitting software
c
c The routine reconstructs the mode shape from the fitted weights for each
c   of the basis functions
      use aarrays
      implicit none
c
c calling argument
      integer imode, j, k, m
      real r( k ), re( k ), im( k )
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
      common / sphlmn / l, mm, n, mode
      integer l, mm, n, mode
c
c externals
      complex sphfn
      real*8 bsjfun, sigraj
c
c local array
      integer ifun( mostn + 1 )
c
c local variables
      complex f
      integer ip, ir
      real ap, cp, p, sp, t, x, y
      real*8 kr, rad
      include 'inc/pi.f'
c
c set scaling
      if( sfpc )then
c basis functions
        if( basset .eq. 'ablj' )then
          maj = m
          kaj = basis
        else if( basset .ne. 'bess' )then
          call crash( 'MSHAPE', 'Unrecognized basis' )
        end if
        ip = 0
        do ir = 1, lastf
          if( msel( ir ) .eq. m )then
            ip = ip + 1
            if( ip .gt. mostn + 1 )call crash( 'MSHAPE',
     +                                      'Function array too small' )
            ifun( ip ) = ir
          end if
        end do
        do ir = j, k
          rad = r( ir )
          re( ir ) = 0
          im( ir ) = 0
          do ip = jp, kp
            if( basset .eq. 'ablj' )then
              naj = nsel( ifun( ip ) )
              ap = sigraj( rad / maxr )
            else
              kr = rad * dble( ip - 1 + minn ) * deltak
              ap = ip * bsjfun( m, kr ) / pi
            end if
            re( ir ) = re( ir ) + ap * alp( ip, imode )
            im( ir ) = im( ir ) + ap * bet( ip, imode )
          end do
        end do
      else if( danl .or. zanl )then
c density or bending analysis
        do ir = j, k
          re( ir ) = alp( ir, imode )
          im( ir ) = bet( ir, imode )
        end do
      else if( lgsp .or. vlgs )then
c logarithmic spirals
        do ir = j, k
          x = 0
          y = 0
          do ip = jp, kp
            p = .25 * real( m * ( ip - ncp / 2 - 1 ) )
            p = p * log( r( ir ) )
            cp = cos( p )
            sp = sin( p )
            x = x + alp( ip, imode ) * cp + bet( ip, imode ) * sp
            y = y - alp( ip, imode ) * sp + bet( ip, imode ) * cp
          end do
          re( ir ) = x / r( ir )**2
          im( ir ) = y / r( ir )**2
        end do
      else if( sphb )then
        call crash( 'MSHAPE', 'Not programmed for new data type' )
c spherical Bessel functions
        mode = imode
        mm = m
        do ir = j, k
          re( ir ) = 0
          im( ir ) = 0
          n = 0
          l = jp / jnmax + m
          do ip = jp, kp
            n = n + 1
            if( n .gt. jnmax )then
              n = 1
              l = l + 2
            end if
            f = sphfn( r( ir ) / rbess, 0., 0. ) / rbess**3
            re( ir ) = re( ir ) + alp( ip, mode ) * real( f ) +
     +                                  bet( ip, mode ) * aimag( f )
            im( ir ) = im( ir ) + alp( ip, mode ) * aimag( f ) -
     +                                   bet( ip, mode ) * real( f )
          end do
        end do
      else
        call crash( 'MSHAPE', 'Unrecognized data type' )
      end if
c convert to t = tme( jt ) only if the mode decays
      if( freq( 2, imode ) .lt. 0. )then
        t = tme( jt ) - tme( kt )
        p = freq( 1, imode ) * t
        cp = cos( p )
        sp = sin( p )
        ap = exp( freq( 2, imode ) * t )
        do ir = j, k
          x = ap * ( re( ir ) * cp - im( ir ) * sp )
          y = ap * ( re( ir ) * sp + im( ir ) * cp )
          re( ir ) = x
          im( ir ) = y
        end do
      end if
      return
      end
