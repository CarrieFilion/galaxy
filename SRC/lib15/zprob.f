      real*8 function zprob( z0 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns a z value randomly chosen from the desired disc model
c
c calling argument
      real*8 z0
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 ranuni, rannor, zmassf
c
c local array
      real*8 w( 17 )
c
c local variables
      integer ifail, ind, ir
      real*8 err, f, tol, x1, x2
c
      zprob = 0
c skip if zero thickness
      if( z0 .gt. 0. )then
        if( iztyp( icmp ) .eq. 2 )then
c Gaussian distribution
          zprob = rannor( 0.d0, z0 )
        else
c numerical solution
          f = ranuni( 0. )
          x2 = 20
          x1 = -x2
          tol = 1.e-6
          ind = 1
          ir = 0
          ifail = 0
c solve for f = zmassd( z / z0 )
          do while ( ind .ne. 0 )
            call fndzro( x1, x2, err, tol, ir, w, ind, ifail )
            err = zmassf( x1 ) - f
          end do
c scale to required thickness
          zprob = z0 * x1
        end if
      end if
      return
      end
