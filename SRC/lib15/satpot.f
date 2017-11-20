      real function satpot( is )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c return the potential arising from an external satellite as the location of
c   of a particle.  Returned value is in internal program units.
c
c calling argument
      integer is
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, n
      real d2
c
c no potential if there is no perturber
      satpot = 0
      if( pertbn )then
c distance of perturber from this particle
        d2 = eptbr * eptbr
        n = 3
        if( twod )n= 2
        do i = 1, n
          d2 = d2 + ( oldc( i, is ) / lscale - xptrb( i ) )**2
        end do
c Plummer sphere potential
        satpot = -mptbr / ( sqrt( d2 ) * gvfac**2 )
      end if
      return
      end
