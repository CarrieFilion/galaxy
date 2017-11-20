      subroutine drmode( mact, nmodes, start, end )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to contour best fitting modes found by modefit
c   This routine is not quite such a mess as it used to be!
c
c calling arguments
      integer mact, nmodes
      real start, end
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      integer nrng
      parameter ( nrng = 501 )
      real rad( nrng )
      common / trcdent / rad
c
c externals
      external cdentr
      real grofu
c
c local arrays
      integer lwrk, nc, nspokes
      parameter ( nspokes = 320, nc = 10 )
      parameter ( lwrk = ( nspokes + 1 ) * nrng )
      real amp( nrng ), cont( nc ), phase( nrng ), work( lwrk )
      real re( nrng ), im( nrng )
      real*8 res( 4 )
c
c local variables
      integer ia, ic, imode, ir, jr, kr, n
      real ampm, ph, r, rmx
      include 'inc/pi.f'
c
      rmx = rgrid( jgrid ) / lscale
      if( p2d .or. p3d )then
c check space
        if( ( nr( jgrid ) .ge. nrng ) .or. ( na .gt. nspokes ) )then
          call jsebuf
          print 200, nrng - 1, nspokes, nr( jgrid ), na
  200     format( ' Insufficient work space - array dimensioned for',
     +            2i4, ' but grid is', 2i4 / ' Error in DRMODE' )
          return
        end if
      end if
c insert arbitrary values
      if( na .le. 0 )na = 80
      if( alpha .le. 0. )alpha = 2 * pi / real( na )
c work over modes
      imode = 0
   10 imode = imode + 1
      if( imode .gt. nmodes )then
        call jschsz( .3 )
        return
      end if
      if( mact .eq. 0 )then
        call rdmode( imode, mact, .false. )
        go to 10
      end if
c draw frame and write header
      call jspage
c      call jschsz( .2 )
      call jssize( .15, .95, .1, .9 )
      call jsescl( -rmx, rmx, -rmx, rmx )
      if( imode .eq. 1 )then
        call jsbldt( 'Run no' )
        call jsbldi( irun, 5 )
        call jsbldt( 'from time' )
        call jsbldf( start, 6, 1 )
        call jsbldt( 'to' )
        call jsbldf( end, 6, 1 )
        if( danl )call jsbldt( ' density' )
        if( zanl )call jsbldt( ' z displacement' )
        if( lgsp )call jsbldt( ' log spi (density)' )
        if( vlgs )then
          if( avlgs )then
            call jsbldt( ' log spi (azimuthal velocity)' )
          else
            call jsbldt( ' log spi (radial velocity)' )
          end if
        end if
        if( sfpc )then
          if( basset .eq. 'ablj' )then
            call jsbldt( 'Abel-Jacobi basis set' )
            call jsbldi( basis, 2 )
          else if( basset .eq. 'bess' )then
            call jsbldt( 'Bessel fn expansion with \\delta_{k}=' )
            call jsbldf( deltak, 4, 2 )
          else
            call crash( 'DRMODE', 'Unrecognized basis' )
          end if
        end if
        call jswrit( -rmx, 1.15 * rmx )
      end if
      if( sphb )then
        if( .not. gronly )then
          call jsbldt( 'patt sp' )
          call jsbldf( freq( 1, imode ), 6, 3 )
        end if
        if( .not. ptonly )then
          call jsbldt( 'grw rt' )
          call jsbldf( freq( 2, imode ), 6, 3 )
        end if
        call jswrit( -rmx, 1.05 * rmx )
        call project( imode, mact )
      else
        call jsaxis( 'x', ' ', 1 )
        call jsaxis( 'y', ' ', 1 )
c draw resonance circles
        if( .not. zanl )then
          call resrad( dble( freq( 1, imode ) ), mact, res )
          do ir = 1, 4
            if( ir .eq. 2 )call jsdash( 0, 1, 0, 1 )
            if( ( res( ir ) .gt. 0. ) .and. ( res( ir ) .lt. rmx ) )then
              call jscirc( 0., 0., sngl( res( ir ) ) )
            end if
          end do
          call jsdash( 0, 0, 0, 0 )
        end if
        if( danl .or. zanl )then
          jr = jp
          kr = kp
          do ir = jr, kr
            if( danl )rad( ir ) = grofu( real( ir - 1 ) ) / lscale
            if( zanl )rad( ir ) =
     +                          ( real( ir ) - .5 ) / ( drfac * lscale )
          end do
c determine interesting radius range
        else
          if( sfpc )then
            radm = maxr / lscale
            rmin = radm / real( nrng )
          else if( lgsp .or. vlgs )then
c radial range for log spirals
            if( res( 3 ) .gt. 0. )then
              rmin = res( 3 )
            else
              rmin = 0.3 * res( 1 )
            end if
            radm = rmx
            if( res( 2 ) .gt. 0.d0 )radm = min( sngl( res( 2 ) ), rmx )
            r = .1 * res( 1 )
            if( r .le. 0. )r = .05 * rtrunc( 1 )
            rmin = min( r, .05 * rtrunc( 1 ) )
          else
            call crash( 'DRMODE', 'Unrecognized data type' )
          end if
          jr = 1
          kr = nrng
          do ir = jr, kr
            rad( ir ) = rmin +
     +                  ( radm - rmin ) * real( ir - 1 ) / real( nrng )
          end do
        end if
c get mode shape
        call mshape( imode, mact, jr, kr, rad, re, im )
c find maximum amplitude
        ampm = 0.
        do ir = jr, kr
          amp( ir ) = sqrt( re( ir )**2 + im( ir )**2 )
          ph = atan2( im( ir ), re( ir ) )
          if( mact .gt. 0 )then
            phase( ir ) = ph / real( mact )
          else
            phase( ir ) = 0
          end if
          ampm = max( amp( ir ), ampm )
        end do
        if( .not. gronly )then
          call jsbldt( 'patt sp' )
          call jsbldf( freq( 1, imode ), 6, 3 )
        end if
        if( .not. ptonly )then
          call jsbldt( 'grw rt' )
          call jsbldf( freq( 2, imode ), 6, 3 )
        end if
        call jsbldt( 'amp scale' )
        call jsbldf( ampm, 8, 4 )
        if( sfpc )then
          call jsbldt( 'maxn' )
          call jsbldi( kp - 1, 2 )
          call jsbldt( 'maxr' )
          call jsbldf( minmaxr / lscale, 10, 2 )
        end if
        call jswrit( -rmx, 1.05 * rmx )
c work over radius and angle
        n = 0
        do ir = jr, kr
          rad( ir - jr + 1 ) = rad( ir )
          do ia = 1, na + 1
            ph = ia - 1
            ph = real( mact ) * ( phase( ir ) - alpha * ph )
            n = n + 1
            work( n ) = amp( ir ) * cos( ph )
          end do
        end do
c contour mode
        do ic = 1, nc
          cont( ic ) = ampm * ( real( ic ) - .5 ) / real( nc )
        end do
        call jsthik( 1 )
        call jscont( work, na + 1, kr - jr + 1, cont, nc, cdentr )
        call jsthik( 3 )
      end if
      go to 10
      end
