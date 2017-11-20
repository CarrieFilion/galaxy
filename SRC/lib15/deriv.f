      real function deriv( x, funct )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c estimates derivative of the given function at the point x
c
c calling arguments
      real x, funct
c
c local variables
      integer nfail, niter
      real anorm, best, dern, dero, diffl, diffn, diffo, d1, d2, d3, d4
      real f0, h
c
c initial step size
      h = 1.
      h = .01 * max( h, abs( x ) )
c decide whether to normalize
      f0 = funct( x )
      anorm = 1
      if( abs( f0 ) .gt. 1. )anorm = 1. / abs( f0 )
c fifth order estimate
      d1 = anorm * ( funct( x + h ) - funct( x - h ) )
      d2 = anorm * ( funct( x + 2. * h ) - funct( x - 2. * h ) )
      d3 = anorm * ( funct( x + 3. * h ) - funct( x - 3. * h ) )
      dero = ( 45. * d1 - 9. * d2 + d3 ) / ( 60. * h )
c seventh order estimate
      d4 = anorm * ( funct( x + 4. * h ) - funct( x - 4. * h ) )
      dern = ( 672. * d1 - 168. * d2 + 32. * d3 - 3. * d4 )
     +                                                    / ( 840. * h )
      best = dern
c monitor convergence
      diffn = abs( dero - dern )
      niter = 1
      nfail = 0
      diffl = diffn
    1 diffo = diffn + 1
c halve step and try again
      do while ( diffn .lt. diffo )
        niter = niter + 1
        dero = dern
        diffo = diffn
        h = .5 * h
        d4 = d2
        d2 = d1
        d3 = anorm * ( funct( x + 3. * h ) - funct( x - 3. * h ) )
        d1 = anorm * ( funct( x + h ) - funct( x - h ) )
        dern = ( 672. * d1 - 168. * d2 + 32. * d3 - 3. * d4 )
     +                                                    / ( 840. * h )
        diffn = abs( dern - dero )
      end do
c check absolute and relative accuracy
      if( ( diffo .gt. 1.e-7 ) .and.
     +    ( diffo .gt. 1.e-5 * abs( dero ) ) )then
c check for reasonable no of iterations
        if( niter .lt. 5 )then
c convergence failed
          nfail = nfail + niter - 1
          niter = 1
c store best result so far
          if( diffo .lt. diffl )then
            diffl = diffo
            best = dero
          end if
          if( nfail .lt. 20 )go to 1
          dero = best
        end if
      end if
c return best estimate
      deriv = dero / anorm
      return
      end

      real*8 function deriv2( x, funct )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c estimates derivative of the given function at the point x
c
c calling arguments
      real*8 x, funct
c
c local variables
      integer nfail, niter
      real*8 anorm, best, dern, dero, diffl, diffn, diffo, d1, d2, d3
      real*8 d4, f0, h
c
c initial step size
      h = 1.
      h = .001 * max( h, abs( x ) )
c decide whether to normalize
      f0 = funct( x )
      anorm = 1
      if( abs( f0 ) .gt. 1.d0 )anorm = 1. / abs( f0 )
c fifth order estimate
      d1 = anorm * ( funct( x + h ) - funct( x - h ) )
      d2 = anorm * ( funct( x + 2. * h ) - funct( x - 2. * h ) )
      d3 = anorm * ( funct( x + 3. * h ) - funct( x - 3. * h ) )
      dero = ( 45. * d1 - 9. * d2 + d3 ) / ( 60. * h )
c seventh order estimate
      d4 = anorm * ( funct( x + 4. * h ) - funct( x - 4. * h ) )
      dern = ( 672. * d1 - 168. * d2 + 32. * d3 - 3. * d4 )
     +                                                    / ( 840. * h )
      best = dern
c monitor convergence
      diffn = abs( dero - dern )
      niter = 1
      nfail = 0
      diffl = diffn
    1 diffo = diffn + 1
c halve step and try again
      do while ( diffn .lt. diffo )
        niter = niter + 1
        dero = dern
        diffo = diffn
        h = .5 * h
        d4 = d2
        d2 = d1
        d3 = anorm * ( funct( x + 3. * h ) - funct( x - 3. * h ) )
        d1 = anorm * ( funct( x + h ) - funct( x - h ) )
        dern = ( 672. * d1 - 168. * d2 + 32. * d3 - 3. * d4 )
     +                                                    / ( 840. * h )
        diffn = abs( dern - dero )
      end do
c check absolute and relative accuracy
      if( ( diffo .gt. 1.d-11 ) .and.
     +    ( diffo .gt. 1.d-13 * abs( dero ) ) )then
c check for reasonable no of iterations
        if( niter .lt. 5 )then
c convergence failed
          nfail = nfail + niter - 1
          niter = 1
c store best result so far
          if( diffo .lt. diffl )then
            diffl = diffo
            best = dero
          end if
          if( nfail .lt. 20 )go to 1
          dero = best
        end if
      end if
c return best estimate
      deriv2 = dero / anorm
      return
      end
