      real*8 function rhoPsi( Psi )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c function needed for Eddington's integral formula DF for an isotropic
c   spherical model.  formula (4-140b) of BT (p 237)
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
c externals
      real*8 phitot, rhohal
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
        err = Psi + phitot( r1 ) - phimax
      end do
      if( ifail .ne. 0 )then
        if( Psi .gt. phitot( 0.d0 ) )then
          r1 = 0
        else
          print *, r1, r2, phimax - phitot( r1 ),
     +             phimax - phitot( r2 ), Psi
          print *, 'IFAIL =', ifail, ' from FNDZRO in RHOPSI'
          call crash( 'RHOPSI', 'Failed to find zero' )
        end if
      end if
c get density
      rhoPsi = rhohal( r1 )
      return
      end
