      real*8 function rcirc( Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c finds the radius of a circular orbit of the given angular momentum
c   in the current potential
c
c calling argument
      real*8 Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 vcirc
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir
      real*8 aLz, res, rl, ru, tol
      include 'inc/pi.f'
c
      rcirc = 0.
      aLz = abs( Lz )
      if( aLz .gt. 0. )then
c set phoney value
        rcirc = -1
c use analytic expression when possible
        if( ncmp .eq. 1 )then
          if( disc( 1 ) )then
c MFK disc - rcirc is a non-monotonic function around outer edge
            if( ( ctype( 1 ) .eq. 'MFK ' ) .and.
     +          ( rmax .le. 1.d0 ) )then
              rcirc = sqrt( aLz / sqrt( 3. * pi / 4. ) )
c MTZ disc
            else if( ctype( 1 ) .eq. 'MTZ ' )then
              rcirc = aLz
c power law disc
            else if( ctype( 1 ) .eq. 'POWE' )then
              res = 1.d0 / ( 1.d0 + dfcns( 3, 1 ) )
              rcirc = aLz**res
            end if
          else
c Jaffe halo
            if( ctype( 1 ) .eq. 'JAFF' )then
              tol = 4. * rscale( icmp ) * cmpmas( icmp ) + aLz**2
              rcirc =
     +             aLz * ( aLz + sqrt( tol ) ) / ( 2. * cmpmas( icmp ) )
            end if
          end if
        end if
      end if
c phoney value unchanged
      if( rcirc .lt. 0.d0 )then
c numerical search for rcirc given by Lz = rcirc * vcirc( rcirc )
        ir = 0
        ind = 1
        ifail = 1
        tol = 1.e-30
c guess range for rcirc
        rl = 0
        ru = max( 10.d0, 2.d0 * rmax )
        ru = max( ru, Lz**4 )
        ru = min( ru, 1.d3 )
        if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'SPLY' ) )ru = 1
        do while ( ind .ne. 0 )
          call fndzro( rl, ru, res, tol, ir, w, ind, ifail )
          res = aLz - rl * vcirc( rl )
        end do
c test for error
        rcirc = rl
        if( ( ifail .ne. 0 ) .and. ( ifail .ne. 5 ) )then
          print *, 'IFAIL = ', ifail, ' from FNDZRO in RCIRC'
          call crash( 'RCIRC', 'Failed to find zero' )
        end if
      end if
      return
      end
