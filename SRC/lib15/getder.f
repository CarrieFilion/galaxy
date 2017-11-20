      subroutine getder( x, absc, ord, n, deriv, ind )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c estimates the first derivative of a function at the point x by reverse
c   communication.  ind must be set to 1 on the first call.  On subsequent
c   calls the value of the ordinate must be set in array ord for each of
c   the n abscissae.  When further subdivision causes no improvement, the
c   best estimate of the derivative is returned with ind set to zero.  The
c   value is generally good to a few parts in 10**4 in single precision
c
c calling arguments
      integer ind, n
      real deriv, ord( 8 ), x, absc( 8 )
c
c local variables
      integer j, k, nfail, niter
      real anorm, best, dern, dero, diffl, diffn, diffo, dx, d1, d2, d3
      real d4, h
      save anorm, best, dern, dero, diffl, diffn, diffo, d1, d2, d3, d4
      save h, nfail, niter
c
c first call - request function value
      if( ind .eq. 1 )then
        n = 1
        absc( 1 ) = x
        ind = 2
c second call - decide whether to normalize and set initial step size
      else if( ind .eq. 2 )then
        anorm = 1
        if( abs( ord( 1 ) ) .gt. 1. )anorm = 1. / abs( ord( 1 ) )
        dx = 1
        if( x .ne. 0. )dx = max( 1., log10( abs( x ) ) )
        h = .01 * dx
c request 8 function values
        n = 8
        do j = 1, 4
          k = 2 * ( j - 1 )
          dx = real( j ) * h
          absc( k + 1 ) = x - dx
          absc( k + 2 ) = x + dx
        end do
        ind = 3
c third call
      else if( ind .eq. 3 )then
c fifth order estimate
        d1 = anorm * ( ord( 2 ) - ord( 1 ) )
        d2 = anorm * ( ord( 4 ) - ord( 3 ) )
        d3 = anorm * ( ord( 6 ) - ord( 5 ) )
        dero = ( 45. * d1 - 9. * d2 + d3 ) / ( 60. * h )
c seventh order estimate
        d4 = anorm * ( ord( 8 ) - ord( 7 ) )
        dern = ( 672. * d1 - 168. * d2 + 32. * d3 - 3. * d4 )
     +                                                    / ( 840. * h )
        best = dern
c monitor convergence
        diffn = abs( dero - dern )
        niter = 1
        nfail = 0
        diffl = diffn
c set arbitrary diffo just to start
    2   diffo = diffn + 1
c halve step and get four new values
    1   h = .5 * h
        n = 4
        do j = 1, 2
          k = 2 * ( j - 1 )
          dx = h
          if( j .eq. 2 )dx = 3 * h
          absc( k + 1 ) = x - dx
          absc( k + 2 ) = x + dx
        end do
        ind = 4
c general call
      else if( ind .eq. 4 )then
c halve step and try again
        niter = niter + 1
        dero = dern
        diffo = diffn
        d4 = d2
        d2 = d1
        d3 = anorm * ( ord( 4 ) - ord( 3 ) )
        d1 = anorm * ( ord( 2 ) - ord( 1 ) )
        dern = ( 672. * d1 - 168. * d2 + 32. * d3 - 3. * d4 )
     +                                                    / ( 840. * h )
        diffn = abs( dern - dero )
c loop back while differences are decreasing
        if( diffn .lt. diffo )go to 1
c check absolute and relative accuracy
        if( ( diffo .gt. 1.e-6 ) .and.
     +      ( diffo .gt. 1.e-4 * abs( dero ) ) )then
c check for a reasonable number of iterations
          if( niter .lt. 5 )then
c convergence failed
            nfail = nfail + niter - 1
            niter = 1
c store best result so far
            if( diffo .lt. diffl )then
              diffl = diffo
              best = dero
            end if
            if( nfail .lt. 20 )go to 2
            dero = best
          end if
        end if
c return best estimate
        deriv = dero / anorm
        ind = 0
      else
        print *, 'ind =', ind
        call crash( 'getder', 'nonsense value of ind input' )
      end if
      return
      end

      subroutine getder2( x, absc, ord, n, deriv, ind )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c estimates the first derivative of a function at the point x by reverse
c   communication.  ind must be set to 1 on the first call.  On subsequent
c   calls the value of the ordinate must be set in array ord for each of
c   the n abscissae.  When further subdivision causes no improvement, the
c   best estimate of the derivative is returned with ind set to zero.  The
c   value is generally good to a few parts in 10**10 in double precision
c
c calling arguments
      integer ind, n
      real*8 deriv, ord( 8 ), x, absc( 8 )
c
c local variables
      integer j, k, nfail, niter
      real*8 anorm, best, dern, dero, diffl, diffn, diffo, dx, d1, d2
      real*8 d3, d4, h
      save anorm, best, dern, dero, diffl, diffn, diffo, d1, d2, d3, d4
      save h, nfail, niter
c
c first call - request function value
      if( ind .eq. 1 )then
        n = 1
        absc( 1 ) = x
        ind = 2
c second call - decide whether to normalize and set initial step size
      else if( ind .eq. 2 )then
        anorm = 1
        if( abs( ord( 1 ) ) .gt. 1.d0 )anorm = 1. / abs( ord( 1 ) )
        dx = 1
        if( x .ne. 0. )dx = max( 1.d0, log10( abs( x ) ) )
        h = .01d0 * dx
c request 8 function values
        n = 8
        do j = 1, 4
          k = 2 * ( j - 1 )
          dx = real( j ) * h
          absc( k + 1 ) = x - dx
          absc( k + 2 ) = x + dx
        end do
        ind = 3
c third call
      else if( ind .eq. 3 )then
c fifth order estimate
        d1 = anorm * ( ord( 2 ) - ord( 1 ) )
        d2 = anorm * ( ord( 4 ) - ord( 3 ) )
        d3 = anorm * ( ord( 6 ) - ord( 5 ) )
        dero = ( 45.d0 * d1 - 9.d0 * d2 + d3 ) / ( 60.d0 * h )
c seventh order estimate
        d4 = anorm * ( ord( 8 ) - ord( 7 ) )
        dern = ( 672.d0 * d1 - 168.d0 * d2 + 32.d0 * d3 - 3.d0 * d4 )
     +                                                  / ( 840.d0 * h )
        best = dern
c monitor convergence
        diffn = abs( dero - dern )
        niter = 1
        nfail = 0
        diffl = diffn
c set arbitrary diffo just to start
    2   diffo = diffn + 1
c halve step and get four new values
    1   h = .5d0 * h
        n = 4
        do j = 1, 2
          k = 2 * ( j - 1 )
          dx = h
          if( j .eq. 2 )dx = 3.d0 * h
          absc( k + 1 ) = x - dx
          absc( k + 2 ) = x + dx
        end do
        ind = 4
c general call
      else if( ind .eq. 4 )then
c halve step and try again
        niter = niter + 1
        dero = dern
        diffo = diffn
        d4 = d2
        d2 = d1
        d3 = anorm * ( ord( 4 ) - ord( 3 ) )
        d1 = anorm * ( ord( 2 ) - ord( 1 ) )
        dern = ( 672.d0 * d1 - 168.d0 * d2 + 32.d0 * d3 - 3.d0 * d4 )
     +                                                  / ( 840.d0 * h )
        diffn = abs( dern - dero )
c loop back while differences are decreasing
        if( diffn .lt. diffo )go to 1
c check absolute and relative accuracy
        if( ( diffo .gt. 1.d-11 ) .and.
     +      ( diffo .gt. 1.d-13 * abs( dero ) ) )then
c check for a reasonable number of iterations
          if( niter .lt. 5 )then
c convergence failed
            nfail = nfail + niter - 1
            niter = 1
c store best result so far
            if( diffo .lt. diffl )then
              diffl = diffo
              best = dero
            end if
            if( nfail .lt. 20 )go to 2
            dero = best
          end if
        end if
c return best estimate
        deriv = dero / anorm
        ind = 0
      else
        print *, 'ind =', ind
        call crash( 'getder', 'nonsense value of ind input' )
      end if
      return
      end
