      real*8 function sderiv2( x, funct, lb, ub )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns an estimate of the gradient of the given real*8 function at
c   the point x.  When x is close to the boundaries, it uses one-sided
c   differnces to avoid seeking values outside the interval lb -> ub
c
c calling arguments
      real*8 lb, ub, x, funct
c
c externals
      external funct
      real*8 deriv2
c
c local variables
      integer nfail, niter
      real*8 anorm, best, dern, dero, diffl, diffn, diffo, f0, f1, f2
      real*8 f3, f4, h
c
      if( ( x .lt. lb ) .or. ( x .gt. ub ) )then
        print *, lb, x, ub
        call crash( 'SDERIV2', 'Argument outside allowed range' )
      end if
c initial step size
      h = 1
      h = .001 * max( h, abs( x ) )
c test whether the argument is far enough from the bounds for central diffs
      diffn = min( x - lb, ub - x )
      if( diffn .gt. h )then
        sderiv2 = deriv2( x, funct )
      else
c determine which bound is closest
        if( x - lb .gt. ub - x )h = -h
c decide whether to normalize
        f0 = funct( x )
        if( abs( f0 ) .gt. 1.d0 )then
          anorm = 1. / f0
          f0 = 1
        else
          anorm = 1
        end if
c second order estimate
        f1 = anorm * funct( x + h )
        f2 = anorm * funct( x + 2. * h )
        dero = .5 * ( -3. * f0 + 4. * f1 - f2 ) / h
c fourth order estimate
        f3 = anorm * funct( x + 3. * h )
        f4 = anorm * funct( x + 4. * h )
        dern = ( -25. * f0 + 48. * f1 - 36. * f2 + 16. * f3 - 3. * f4 )
     +         / ( 12. * h )
        best = dern
c monitor convergence
        diffn = abs( dero - dern )
        niter = 1
        nfail = 0
        diffl = diffn
    1   diffo = 1.01 * diffn
c halve step and try again
        do while ( diffn .lt. diffo )
          niter = niter + 1
          dero = dern
          diffo = diffn
          h = .5 * h
          f4 = f2
          f2 = f1
          f1 = anorm * funct( x + h )
          f3 = anorm * funct( x + 3. * h )
         dern = ( -25. * f0 + 48. * f1 - 36. * f2 + 16. * f3 - 3. * f4 )
     +           / ( 12. * h )
          diffn = abs( dern - dero )
        end do
c absolute or relative accuracy
        if( ( diffo .gt. 1.d-10 ) .and.
     +      ( diffo .gt. 1.d-12 * abs( dero ) ) )then
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
            if( nfail .lt. 10 )go to 1
            dero = best
          end if
        end if
c return best estimate
        sderiv2 = dero / anorm
      end if
      return
      end
