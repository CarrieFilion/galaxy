      integer function zoneno( is, new )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c determines current zone number of particle based on its newc position
c
c calling argument
      integer is
      logical new
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local variables
      integer i
      real rad, x, y, z
c
c skip if nothing to do
      zoneno = 1
      if( ( nzones .gt. 1 ) .or. ( nguard .gt. 0 ) )then
c find radius of particle
        if( new )then
c          zoneno = max( iz( is ), 1 )
c          zoneno = min( zoneno, nzones )
c should really use xcpred, but that supposes the zone is already known!
          x = newc( 1, is ) - xcen( 1, jgrid )
          y = newc( 2, is ) - xcen( 2, jgrid )
          if( threed )z = newc( 3, is ) - xcen( 3, jgrid )
        else
          x = oldc( 1, is ) - xcen( 1, jgrid )
          y = oldc( 2, is ) - xcen( 2, jgrid )
          if( threed )z = oldc( 3, is ) - xcen( 3, jgrid )
        end if
        if( sphbnd )then
          rad = sqrt( x * x + y * y + z * z )
        else
          rad = sqrt( x * x + y * y )
        end if
c find current zone of particle
        if( nzones .gt. 1 )then
          do i = 1, nzones
            if( rad .ge. zbr( i ) )zoneno = i
          end do
        end if
        if( ( zoneno .eq. 1 ) .and. ( nguard .gt. 0 ) )then
          do i = 1, nguard
            if( rad .lt. gbr( i ) )zoneno = 1 - i
          end do
        end if
      end if
      return
      end
