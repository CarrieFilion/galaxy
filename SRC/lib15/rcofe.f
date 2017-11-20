      real*8 function rcofe( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns radius of circular orbit with energy E
c
c calling argument
      real*8 E
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 phitot, vcirc
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir
      real*8 a, b, c, cc, err, rout, tol, xi
      include 'inc/pi.f'
c
      rcofe = -1
c check for particle with minimum allowed energy
      if( E .le. Emine )then
        rcofe = rhole
c check whether solution exists
        a = E - Emine
        if( Emine .ne. 0 )a = a / Emine
        if( abs( a ) .gt. 1.d-5 )then
          print *, E, Emine, E - Emine
          call crash( 'RCOFE', 'Unattainable energy' )
        end if
      end if
c use analytic expression when possible
      if( ( rcofe .lt. 0.d0 ) .and. ( ncmp .eq. 1 ) )then
        if( disc( 1 ) )then
          cc = dfcns( 3, 1 )
c Maclaurin/Freeman/Kalnajs disc
          if( ctype( 1 ) .eq. 'MFK ' )then
            rcofe = sqrt( 1. + 4. * E / ( 3. * pi ) )
c MTZ disc - circular orbit given by E = ln( r ) + .5
          else if( ctype( 1 ) .eq. 'MTZ ' )then
            rcofe = exp( E - .5 )
c power law disc - E = r^{2c}/(2c) + r^{2c}/2
          else if( ctype( 1 ) .eq. 'POWE' )then
            tol = .5d0 / cc
            rcofe = ( 2.d0 * E / ( 1.d0 + 1.d0 / cc ) )**tol
c Sawamura's polynomial discs
          else if( ctype( 1 ) .eq. 'SPLY' )then
            if( cc .eq. 0. )then
              rcofe = sqrt( 1. + 4. * E / ( 3. * pi ) )
            else
              a = 3. * cc
              b = 2. * ( 1. - cc )
              c = cc / 3. + E * ( 1. + .4 * cc ) / ( .375 * pi )
              xi = ( -b + sqrt( b * b - 4. * a * c ) ) / ( 2. * a )
              rcofe = sqrt( max( 1. - xi, 0.d0 ) )
            end if
          end if
        else
          if( ctype( 1 ) .eq. 'KEPL' )then
c Keplerian rotation curve
            rcofe = -cmpmas( 1 ) / ( 2. * E )
          else if( ctype( 1 ) .eq. 'UNIS' )then
c uniform sphere
            a = 2. * E * rscale( 1 ) / cmpmas( 1 )
            if( a .le. -1.d0 )then
              rcofe = rscale( 1 ) * sqrt( .5 * ( 3. + a ) )
            else
              rcofe = -cmpmas( 1 ) / ( 2. * E )
            end if
          end if
        end if
      end if
c none of the above
      if( rcofe .lt. 0.d0 )then
c guess initial range
        rcofe = 0.
        rout = 2
        ifail = 1
c expand initial range until a zero is found
        do while ( ifail .eq. 1 )
          rout = 10. * rout
          tol = 1.e-30
          ir = 1
          ind = 1
          ifail = 1
c find zero within current range
          do while ( ind .ne. 0 )
            call fndzro( rcofe, rout, err, tol, ir, w, ind, ifail )
            err = E - ( phitot( rcofe ) + .5 * vcirc( rcofe )**2 )
          end do
        end do
        if( ( ifail .ne. 0 ) .and. ( ifail .lt. 4 ) )then
          print *, 'IFAIL =', ifail, ' from FNDZRO in RCOFE'
          call crash( 'RCOFE', 'Failed to find zero' )
        end if
      end if
      return
      end
