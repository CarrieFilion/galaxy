      real*8 function rofPsi( Psi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c function needed for Eddington inversion
c
c calling argument
      real*8 Psi
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 Phitot
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind, ir
      real*8 err, r1, r2, tol
      include 'inc/pi.f'
c
c guess initial range
      r1 = 0.
      r2 = 1.5 * rmax
      tol = 1.d-10
      ir = 1
      ind = 1
      ifail = 1
c find radius for this Psi
      do while ( ind .ne. 0 )
        call fndzro( r1, r2, err, tol, ir, w, ind, ifail )
        err = Psi + Phitot( r1 ) - Phimax
      end do
      if( ifail .ne. 0 )then
        if( Psi .gt. Phitot( 0.d0 ) )then
          r1 = 0
        else
          print *, r1, r2, Phimax - Phitot( r1 ),
     +             Phimax - Phitot( r2 ), Psi
          print *, 'IFAIL =', ifail, ' from FNDZRO in RHOPSI'
          call crash( 'ROFPSI', 'Failed to find zero' )
        end if
      end if
      rofPsi = r1
      return
      end
