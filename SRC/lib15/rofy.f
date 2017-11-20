      real*8 function rofy( y )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Evaluates r for a given value of y.  y = r * phi( r ) / ( rstar * phi( 0 ) )
c   y is used in definition of Kalnajs's distribution functions
c   Ap J 205, p751 (1976)
c
c calling argument
      real*8 y
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 phitot
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir
      real*8 r1, tol, yofr
c
      if( ( ncmp .gt. 1 ) .or. ( .not. disc( 1 ) )
     +                          )call crash( 'ROFY', 'zero disc model' )
c check for y in valid range
      rofy = -1
      if( ( y .ge. 0. ) .and. ( y .lt. 1. ) )then
c Kuz'min/Toomre disc
        if( ( ctype( 1 ) .eq. 'KT  ' ) .or.
     +      ( ctype( 1 ) .eq. 'SOFT' ) )then
          rofy = 1. - y * y
          rofy = rstar * y / sqrt( rofy )
c Isochrone disc
        else if( ctype( 1 ) .eq. 'ISOC' )then
          rofy = 2. * y / ( 1. - y * y )
        end if
      end if
c general expression
      if( rofy .lt. 0.d0 )then
c numerical search for zero of y( r ) - y
        tol = 1.e-8
        ir = 0
        ind = 1
        ifail = 1
c set range for search
        rofy = 0.
        r1 = 1000.
c find zero
        do while ( ind .ne. 0 )
          call fndzro( rofy, r1, yofr, tol, ir, w, ind, ifail )
          yofr = rofy * phitot( rofy ) / ( phi0 * rstar ) - y
        end do
c check for errors
        if( ifail .gt. 0 )then
          if( ifail .gt. 1 )then
            print *, 'ifail = ', ifail, ' from FNDZRO for arg ', y
            call crash( 'ROFY', 'Failed to find zero' )
          end if
c check for valid calling argument
          if( ( y .lt. -.02 ) .or. ( y .gt. 1.02 ) )then
            print *, 'Argument was ', y
            call crash( 'ROFY', 'Argument out of valid range' )
          end if
c fudge for small extension to range of y (needed for derivatives)
          if( y .lt. 0. )then
            rofy = 0.
          else if( y .ge. 1. )then
            rofy = 10.**( 32. / dfcns( 1, 1 ) )
          else
c r is very large - use rough estimate
            rofy = 1. / sqrt( 1. - y )
          end if
        end if
      end if
      return
      end
