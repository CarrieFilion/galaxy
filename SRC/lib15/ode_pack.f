      subroutine ode_int( x1, x2, nvar, ystart, eps, derivs, ifail )
      implicit none
c adapted from Numerical Recipes p559 to replace NAG d02baf
c
c calling arguments
c         x1, x2         integration range
c         nvar           number of ODEs
c         ystart( nvar ) initial y values
c         eps            requested accuracy
c         derivs         subroutine that returns yprime at coordinate x
c         ifail          set to zero for hard fail, 1 for soft fail
c
      integer ifail, nvar
      external derivs
      real*8 eps, x1, x2, ystart( nvar )
c
c local arrays
      real*8, allocatable :: dydx( : ), y( : ), yscal( : )
c
c local variables
      integer i, maxstp, nbad, nok, nstp
      real*8 h, hdid, hnext, tiny, x
      parameter ( maxstp = 10000, tiny = 1.d-30 )
c
c check calling arguments
      if( abs( x2 - x1 ) .lt. 1.d-30 )then
        if( ifail .eq. 0 )call crash( 'ODE_INT', 'Range too small' )
        ifail = 1
        return
      end if
c
      allocate ( dydx( nvar ), y( nvar ), yscal( nvar ) )
c
      x = x1
c guess initial step
      h = .0005 * ( x2 - x1 )
      nok = 0
      nbad = 0
      do i = 1, nvar
        y( i ) = ystart( i )
      end do
c main step loop
      do nstp = 1, maxstp
        call derivs( x, y, dydx )
c scaling to monitor accuracy - can be modified if need be
        do i = 1, nvar
          yscal( i ) = abs( y( i ) ) + abs( h * dydx( i ) ) + tiny
        end do
c limit last step
        if( ( x + h - x2 ) * ( x + h - x1 ) .gt. 0. )h = x2 - x
c advance
        call nrrkqs(
     +        y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs )
        if( hdid .ne. h )then
          nbad = nbad + 1
        else
          nok = nok + 1
c test whether done
          if( ( x - x2 ) * ( x2 - x1 ) .ge. 0. )then
            do i = 1, nvar
              ystart( i ) = y( i )
            end do
            return
          end if
        end if
c step size for next step
        h = hnext
      end do
      if( ifail .eq. 0 )then
        print *, 'nbad, nok', nbad, nok
        call crash( 'ODE_INT', 'too many steps needed' )
      end if
      ifail = 2
      return
      end

      subroutine ode_toy( x1, x2, nvar, ystart, eps, xint, m, yt,
     +                    derivs, ifail )
      implicit none
c adapted from Numerical Recipes p559 to replace NAG d02bgf
c
c calling arguments
c         x1, x2         integration range
c         nvar           number of ODEs
c         ystart( nvar ) initial y values
c         eps            requested accuracy
c         xint           ignored here - included for consistency with NAG
c         m              the y coordinate to be monitored
c         yt             the value requested
c         derivs         subroutine that returns yprime at coordinate x
c         ifail          set to zero for hard fail, 1 for soft fail
c
      integer ifail, m, nvar
      external derivs
      real*8 eps, xint, x1, x2, ystart( nvar ), yt
c
c external
      real*8 aitken2
c
c local arrays
      real*8, allocatable :: dydx( : ), y( : ), yscal( : )
c
c local variables
      integer i, iextra, maxstp, nbad, nok, ns, nstp
      logical done
      real*8 h, hdid, hnext, tiny, x, xl( 6 ), yl( 6 )
      parameter ( maxstp = 10000, tiny = 1.d-30 )
c
      allocate ( dydx( nvar ), y( nvar ), yscal( nvar ) )
c
c check calling arguments
      if( m .gt. nvar .or. m .lt. 1 )then
        if( ifail .eq. 0 )call crash( 'ODE_TOY', 'Invalid m' )
        ifail = 1
        return
      end if
      if( abs( x2 - x1 ) .lt. 1.d-20 )then
        if( ifail .eq. 0 )call crash( 'ODE_TOY', 'Range too small' )
        ifail = 1
        return
      end if
c
      x = x1
c guess initial step
      h = .0005 * ( x2 - x1 )
      nok = 0
      nbad = 0
      xl( 1 ) = x1
      yl( 1 ) = ystart( m )
      ns = 1
      done = .false.
      iextra = 0
      do i = 1, nvar
        y( i ) = ystart( i )
      end do
c main step loop
      do nstp = 1, maxstp
        call derivs( x, y, dydx )
c scaling to monitor accuracy - can be modified if need be
        do i = 1, nvar
          yscal( i ) = abs( y( i ) ) + abs( h * dydx( i ) ) + tiny
        end do
c advance
        call nrrkqs(
     +        y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs )
        if( hdid .ne. h )then
          nbad = nbad + 1
        else
          nok = nok + 1
c test whether done
          if( .not. done )done =
     +                   ( ( y( m ) - yt ) * ( yl( ns ) - yt ) .le. 0. )
          if( iextra .eq. 2 )then
c return value of x where y( m ) = yt
            ns = ns + 1
            xl( ns ) = x
            yl( ns ) = y( m )
c interpolation for unequally spaced abscissae
            x1 = aitken2( yl, xl, ns, yt )
            ifail = 0
            return
          end if
c need 2 extra steps to improve interpolation
          if( done )iextra = iextra + 1
c keep last 5 x,y values in a table so result can be obtained by interpolation
          if( ns .lt. 5 )then
            ns = ns + 1
            xl( ns ) = x
            yl( ns ) = y( m )
          else
            do i = 1, ns - 1
              xl( i ) = xl( i + 1 )
              yl( i ) = yl( i + 1 )
            end do
            xl( ns ) = x
            yl( ns ) = y( m )
          end if
        end if
c step size for next step
        h = hnext
      end do
      if( ifail .eq. 0 )then
        print *, 'nbad, nok', nbad, nok
        call crash( 'ODE_TOY', 'too many steps needed' )
      end if
      ifail = 2
      return
      end

      subroutine ode_tab( x1, xend, nvar, ystart, derivs, eps, output,
     +                    ifail )
      implicit none
c adapted from Numerical Recipes p559 to replace NAG d02cjf
c
c calling arguments
c         x1, xend       integration range
c         nvar           number of ODEs
c         ystart( nvar ) initial y values
c         derivs         subroutine that returns yprime at coordinate x
c         eps            requested accuracy
c         output         an external subroutine to store intermediate results
c         ifail          set to zero for hard fail, 1 for soft fail
c
      integer ifail, nvar
      external derivs, output
      real*8 eps, x1, xend, ystart( nvar )
c
c local arrays
      real*8, allocatable :: dydx( : ), y( : ), yscal( : )
c
c local variables
      integer i, maxstp, nbad, nok, nstp
      real*8 h, hdid, hnext, tiny, x, x2
      parameter ( maxstp = 10000, tiny = 1.d-30 )
c
      allocate ( dydx( nvar ), y( nvar ), yscal( nvar ) )
c
c check calling arguments
      if( abs( xend - x1 ) .lt. 1.d-20 )then
        if( ifail .eq. 0 )call crash( 'ODE_TAB', 'Range too small' )
        ifail = 1
        return
      end if
c
      x = x1
c guess initial step
      h = .0005 * ( xend - x1 )
      nok = 0
      nbad = 0
      do i = 1, nvar
        y( i ) = ystart( i )
      end do
c get first stop point
      x2 = x1
      call output( x2, y )
      if( ( x2 - x1 ) .lt. 1.d-20 )then
        if( ifail .eq. 0
     +       )call crash( 'ODE_TAB', 'First stop point already passed' )
        ifail = 3
        return
      end if
      if( h .gt. 0.d0 )then
        x2 = min( x2, xend )
      else
        x2 = max( x2, xend )
      end if
c main step loop
      do nstp = 1, maxstp
        call derivs( x, y, dydx )
c scaling to monitor accuracy - can be modified if need be
        do i = 1, nvar
          yscal( i ) = abs( y( i ) ) + abs( h * dydx( i ) ) + tiny
        end do
c limit step size to end at x2
        if( ( x + h - x2 ) * ( x + h - x1 ) .gt. 0. )h = x2 - x
c advance
        call nrrkqs(
     +        y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs )
        if( hdid .ne. h )then
          nbad = nbad + 1
        else
          nok = nok + 1
c test whether done
          if( ( x - x2 ) * ( x2 - x1 ) .ge. 0. )then
c save intermediate result and get next stop point
            call output( x2, y )
c test whether done
            if( ( h .gt. 0.d0 .and. x2 .gt. xend ) .or.
     +          ( h .lt. 0.d0 .and. x2 .lt. xend ) )then
              do i = 1, nvar
                ystart( i ) = y( i )
              end do
              ifail = 0
              return
            end if
          end if
        end if
c step size for next step
        h = hnext
      end do
      if( ifail .eq. 0 )then
        print *, 'nbad, nok', nbad, nok
        call crash( 'ODE_TAB', 'too many steps needed' )
      end if
      ifail = 2
      return
      end

      subroutine nrrkqs( y, dydx, n, x, htry, eps, yscal, hdid, hnext,
     +                   derivs )
c adapted from Numerical Recipes p558
      implicit none
c
c calling arguments
      integer n
      external derivs
      real*8 eps, hdid, hnext, htry, x, dydx( n ), y( n ), yscal( n )
c
c local arrays
      real*8, allocatable, save :: yerr( : ), ytemp( : )
c
c local variables
      integer i, isave
      real*8 errmax, h, xnew, safety, pgrow, pshrnk, errcon
      parameter ( safety = 0.9, pgrow = -.2,
     +            pshrnk = -.25, errcon = 1.89d-10 )
c
      data isave / 0 /
c
      if( isave .ne. n ) then
        if( allocated ( yerr ) )deallocate ( yerr, ytemp )
        allocate ( yerr( n ), ytemp( n ) )
        isave = n
      end if
c
      h = htry
    1 call nrrkck( y, dydx, n, x, h, ytemp, yerr, derivs )
      errmax = 0.
      do i = 1, n
        errmax = max( errmax, abs( yerr( i ) / yscal( i ) ) )
      end do
      errmax = errmax / eps
      if( errmax .gt. 1. )then
        h = safety * h * ( errmax**pshrnk )
        if( h .lt. 0.1 * h )then
          h = .1 * h
        endif
        xnew = x + h
        if( xnew .eq. x )call crash( 'NRRKCK', 'stepsize underflow' )
        go to 1
      else
        if( errmax .gt. errcon )then
          hnext = safety * h * ( errmax**pgrow )
        else
          hnext = 5. * h
        end if
        hdid = h
        x = x + h
        do i = 1, n
          y( i ) = ytemp( i )
        end do
        return
      end if
      end

      subroutine nrrkck( y, dydx, n, x, h, yout, yerr, derivs )
c adapted from Numerical Recipes
      implicit none
c
c calling arguments
      integer n
      external derivs
      real*8 h, x, dydx( n ), y( n ), yerr( n ), yout( n )
c
c local arrays
      real*8, allocatable, save ::
     +                  ak2(:), ak3(:), ak4(:), ak5(:), ak6(:), ytemp(:)
c
c local variables
      integer i, isave
      real*8 a2, a3, a4, a5, a6, b21, b31, b32, b41, b42, b43, b51, b52
      real*8 b53, b54, b61, b62, b63, b64, b65, c1, c3, c4, c6, dc1, dc3
      real*8 dc4, dc5, dc6
      parameter ( a2 = .2d0,
     +            a3 = .3d0,
     +            a4 = .6d0,
     +            a5 = 1.d0,
     +            a6 = .875d0,
     +            b21 = .2d0,
     +            b31 = 3.d0 / 40.d0,
     +            b32 = 9.d0 / 40.d0,
     +            b41 = .3d0,
     +            b42 = -.9d0,
     +            b43 = 1.2d0,
     +            b51 = -11.d0 / 54.d0,
     +            b52 = 2.5d0,
     +            b53 = -70.d0 / 27.d0,
     +            b54 = 35.d0 / 27.d0,
     +            b61 = 1631.d0 / 55296.d0,
     +            b62 = 175.d0 / 512.d0,
     +            b63 = 575.d0 / 13824.d0,
     +            b64 = 44275.d0 / 110592.d0,
     +            b65 = 253.d0 / 4096.d0,
     +            c1 = 37.d0 / 378.d0,
     +            c3 = 250.d0 / 621.d0,
     +            c4 = 125.d0 / 594.d0,
     +            c6 = 512.d0 / 1771.d0,
     +            dc1 = c1 - 2825.d0 / 27648.d0,
     +            dc3 = c3 - 18575.d0 / 48384.d0,
     +            dc4 = c4 - 13525.d0 / 55296.d0,
     +            dc5 = -277.d0 / 14336.d0,
     +            dc6 = c6 - .25d0 )
      data isave / 0 /
c
      if( isave .ne. n ) then
        if( allocated ( ak2 )
     +                    )deallocate ( ak2, ak3, ak4, ak5, ak6, ytemp )
        allocate ( ak2( n ), ak3( n ), ak4( n ), ak5( n ),
     +             ak6( n ), ytemp( n ) )
        isave = n
      end if
c
      do i = 1, n
        ytemp( i ) = y( i ) + b21 * h * dydx( i )
      end do
      call derivs( x + a2 * h, ytemp, ak2 )
      do i = 1, n
        ytemp( i ) = y( i ) + h * ( b31 * dydx( i ) + b32 * ak2( i ) )
      end do
      call derivs( x + a3 * h, ytemp, ak3 )
      do i = 1, n
        ytemp( i ) = y( i ) + h * ( b41 * dydx( i ) + b42 * ak2( i ) +
     +               b43 * ak3( i ) )
      end do
      call derivs( x + a4 * h, ytemp, ak4 )
      do i = 1, n
        ytemp( i ) = y( i ) + h * ( b51 * dydx( i ) + b52 * ak2( i ) +
     +               b53 * ak3( i ) + b54 * ak4( i ) )
      end do
      call derivs( x + a5 * h, ytemp, ak5 )
      do i = 1, n
        ytemp( i ) = y( i ) + h * ( b61 * dydx( i ) + b62 * ak2( i ) +
     +               b63 * ak3( i ) + b64 * ak4( i ) + b65 * ak5( i ) )
      end do
      call derivs( x + a6 * h, ytemp, ak6 )
      do i = 1, n
        yout( i ) = y( i ) + h * ( c1 * dydx( i ) + c3 * ak3( i ) +
     +              c4 * ak4( i ) + c6 * ak6( i ) )
      end do
      do i = 1, n
        yerr( i ) = h * ( dc1 * dydx( i ) + dc3 * ak3( i ) +
     +                dc4 * ak4( i ) + dc5 * ak5( i ) + dc6 * ak6( i ) )
      end do
      return
      end
