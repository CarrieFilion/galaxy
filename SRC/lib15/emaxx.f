      real*8 function emaxx( x )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns value of e for a circular orbit in hard KT potential
c
c calling argument
      real*8 x
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      integer ifail, ind, ir
      real*8 a1, a2, a3, f, r1, r2, tol, x2, w( 4 )
c
      emaxx = 0.
      if( ctype( icmp ) .eq. 'KT  ' )then
c Kuz'min/Toomre disc
        if( x .lt. 1.d0 )then
          x2 = x * x
          a1 = 1. - x2
          a2 = 2. - 3. * x2
          a3 = -3. * x2
          ir = 0
          ind = 1
          ifail = 0
          tol = 1.d-6
          r1 = 0
          r2 = 1
          if( x .gt. .9d0 )r2 = 1. / ( 1. - x )
          do while ( ind .ne. 0 )
            call fndzro( r2, r1, f, tol, ir, w, ind, ifail )
            f = a1 * r2**3 + a2 * r2 * r2 + a3 * r2 - x2
          end do
          emaxx = 1. / sqrt( 1. + r2 ) - .5 * r2 / ( 1. + r2 )**1.5
        end if
      else
        call crash( 'EMAXX', 'Unknown disc' )
      end if
      return
      end
