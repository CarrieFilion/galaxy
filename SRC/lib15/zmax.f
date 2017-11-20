      real*8 function zmax( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns maximum vertical excursion at a given r for a particle of
c   fixed E & h
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/orbval.f'
c
c external
      real*8 phisph
c
c local array
      real*8 w( 4 )
c
c local variables
      integer ifail, ind
      real*8 err, target, tol, z1, z2
c
c guess initial range
      z1 = 0
      z2 = 100
c set target value for potential
      target = crrntE
      if( r .gt. 0. )target = target - .5 * ( crrntL / r )**2
c
      ifail = 1
      ind = 1
      tol = 1.d-10
      do while ( ind .ne. 0 )
        call fndzro( z1, z2, err, tol, 0, w, ind, ifail )
        err = target - phisph( r, z1 )
      end do
c return result
      zmax = z1
      if( ifail .ne. 0 )then
        z1 = min( z1, z2 )
c abort only if error is significant
        if( ifail .ne. 4 )then
          err = target - phisph( r, z1 )
          if( abs( err ) .gt. 1.d-8 )then
            print *, 'IFAIL =', ifail, ' from FNDZRO in ZMAX'
            call crash( 'ZMAX', 'could not find zero' )
          end if
c return default value
          zmax = 0
        end if
      end if
      return
      end
