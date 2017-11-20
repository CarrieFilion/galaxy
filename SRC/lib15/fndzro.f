      subroutine fndzro( a, b, resid, tol, ir, w, ind, ifail )
c routine to find a zero of a function by reverse communication
c   it has the same calling arguments and functionality as C05AZF from NAg
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      integer ifail, ind, ir
      real*8 a, b, resid, tol, w( 4 )
c
c externals
      real*8  x02ajf, x02akf
      integer p01abf
      external x02ajf, x02akf, p01abf
c
c local variables
      integer i
      real*8 d, y, z
      save z
c
      i = 0
c first call with no values set
      if( ind .eq. 1 )then
        ind = 2
c second call, or first call with w( 1 ) set to f( b )
      else if( ( ind .eq. 2 ) .or. ( ind .eq. -1 ) )then
c check arguments this one time
        if( ( tol .le. 0.d0 ) .or.
     +      ( ir .lt. 0 ) .or. ( ir .gt. 2 ) )i = 3
        if( ind .eq. 2 )then
          w( 1 ) = resid
          w( 3 ) = a
          w( 4 ) = b
          a = b
        else
c w( 1 ) is set to f( b ), so values are reversed
          w( 4 ) = a
          w( 3 ) = b
        end if
c set values for w( 2 ) and z
        w( 2 ) = -w( 1 )
        z = abs( w( 1 ) )
        ind = 3
      else if( ind .eq. 3 )then
        w( 2 ) = resid
c check that residuals have opposite signs initially
        if( sign( 1.d0, w( 1 ) ) .eq. sign( 1.d0, w( 2 ) ) )then
c set fail flag only if initial residuals are both non-zero
          if( ( w( 1 ) .ne. 0.d0 ) .and. ( w( 2 ) .ne. 0.d0 ) )then
            i = 1
            a = w( 3 )
            b = w( 4 )
          end if
        end if
        ind = 4
      else if( ind .ne. 4 )then
c inappropriate ind
        i = 2
      end if
c generic call
      if( ( ind .eq. 4 ) .and. ( i .eq. 0 ) )then
c finished if either residual is zero
        if( w( 1 ) .eq. 0.d0 )then
          a = w( 3 )
          ind = 0
          resid = 0
          i = 0
        else if( w( 2 ) .eq. 0.d0 )then
          a = w( 4 )
          ind = 0
          resid = 0
          i = 0
        else
c determine which pair of abscissae straddle the zero
          if( sign( 1.d0, w( 1 ) ) .ne. sign( 1.d0, resid ) )then
            w( 2 ) = resid
            w( 4 ) = a
          else
            w( 1 ) = resid
            w( 3 ) = a
          end if
c convergence test
          z = min( z, abs( resid ) )
          y = z
          if( ir .eq. 0 )y = max( 1.d0, z )
          if( ir .eq. 1 )y = 1
          d = abs( w( 3 ) - w( 4 ) )
          if( d .lt. 2. * tol * y )then
            ind = 0
c if last residuals are not the smallest, this could be a pole rather than a 0
            if( min( abs( w( 1 ) ), abs( w( 2 ) ) ) .gt. z )i = 4
c determine closest value
            if( abs( w( 1 ) ) .lt. abs( w( 2 ) ) )then
              a = w( 3 )
              resid = w( 1 )
            else
              a = w( 4 )
              resid = w( 2 )
            end if
          else
c not converged yet - bisect range and try again
            a = .5 * ( w( 3 ) + w( 4 ) )
            if( ( a .eq. w( 3 ) ) .or. ( a .eq. w( 4 ) ) )then
c no further subdivision possible
              ind = 0
              i = 5
            end if
          end if
        end if
      end if
c return flag
      if( i .ne. 0 )then
        if( ifail .eq. 0 )then
          print *, 'IFAIL =', i
          call crash( 'FNDZRO', 'Hard fail requested' )
        else
          ind = 0
        end if
      end if
      if( ind .eq. 0 )ifail = i
      return
      end
