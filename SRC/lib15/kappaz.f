      real function kappaz( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c compute the vertical oscillation frequency in the mid-plane at the requested
c   radius
c both the calling argument and the result are in natural units
c
c calling argument
      real r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local arrays
      real z( 8 ), fz( 8 )
c
c local variables
      integer i, ind, iuse, n
      real rad, fr, zero
      parameter ( zero = 0. )
      save iuse
      data iuse / -1 /
c
      rad = r * lscale
c estimate first derivatve of the vertical force by reverse communication
      ind = 1
      do while ( ind .ne. 0 )
        call getder( zero, z, fz, n, kappaz, ind )
        if( ind .gt. 0 )then
          do i = 1, n
            call mnfrz( rad, z( i ), fr, fz( i ) )
          end do
        end if
      end do
c
      if( kappaz .le. 0. )then
        kappaz = sqrt( -kappaz ) / ts
      else
        if( iuse .ne. istep )then
          print *, 'kappa z imaginary for r =', r, ' at step', istep
          iuse = istep
        end if
        kappaz = 0
      end if
      return
      end
