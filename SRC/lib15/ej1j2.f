      real*8 function Ej1j2( aj1, aj2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the energy of a particle on a planar orbit with the given actions
c
c calling arguments
      real*8 aj1, aj2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 actj1, Emax, Emin
c
c local arrays
      real*8 c( 4 )
c
c local variables
      integer ifail, ind, ir
      logical set
      real*8 err, E2, tol
c
      set = .false.
c use analytic expression when available
      if( ncmp .eq. 1 )then
        if( ctype( 1 ) .eq. 'ISOC' )then
c isochrone disc
          Ej1j2 = -.5 / ( aj1 + .5 * abs( aj2 )
     +                    + sqrt( .25 * aj2 * aj2 + 1 ) )**2
          set = .true.
        end if
      else if( aj1 .le. 0.d0 )then
c circular orbit
        Ej1j2 = Emin( aj2 )
        set = .true.
      end if
      if( .not. set )then
c close to maximum allowed energy
        E2 = Emax( aj2 )
        err = aj1 - actj1( E2, aj2 )
        if( err .gt. 1.d-6  )then
          print *, aj1, actj1( E2, aj2 )
          call crash( 'EJ1J2', 'Impossible calling arguments' )
        else if( err .gt. -1.d-10 )then
          Ej1j2 = E2
c general cases
        else
c set bounds on range of E
          Ej1j2 = Emin( aj2 )
c find zero
          tol = 1.e-30
          ir = 1
          ind = 1
          ifail = 1
          do while ( ind .ne. 0 )
            call fndzro( Ej1j2, E2, err, tol, ir, c, ind, ifail )
            err = aj1 - actj1( Ej1j2, aj2 )
          end do
c check for error
          if( ( ifail .eq. 1 ) .or. ( ifail .eq. 2 ) )then
            print *, 'ifail = ', ifail, ' from FNDZRO in EJ1J2'
            call crash( 'EJ1J2', 'Failed to find zero' )
          end if
        end if
      end if
      return
      end
