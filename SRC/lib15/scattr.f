      subroutine scattr( jst )
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to scatter the current group of particles in / buffer / back
c   to the appropriate locations in / ptcls / and to make links for the
c   current list to which each belongs.
c See the comments in gather for the meaning of the / ptcls / variables
      use aarrays
      implicit none
c
c calling arguments
      integer jst
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
      integer i, ia, is, j, n
      real a
      equivalence ( ia, a )
c
      n = 1
      if( parallel )n = myid + 1
c scatter particles back to storage area
      do is = 1, jst
        j = loc( is )
        do i = 1, ncoor
          ptcls( j + i ) = newc( i, is )
        end do
c particle weight, if not uniform
        if( uqmass )ptcls( j + ncoor + 1 ) = pwt( is )
c population flag
        ia = iflag( is )
        ptcls( j + nwpp - 1 ) = a
c make links for new list - particles off the "grid"
        if( iz( is ) .eq. nlists )then
          i = nlists
        else
c list number determined by current zone and/or plane
          if( nplanes .eq. 1 )then
            i = iz( is ) + nguard
            i = i + ( label( is ) - 1 ) * ( nzones + nguard )
          else
c no zones if multiple planes - no good reason now?
            i = ipln( is )
            i = i + nplanes * ( label( is ) - 1 )
          end if
        end if
        ia = islist( 2, i, n )
        ptcls( j + nwpp ) = a
        islist( 2, i, n ) = j
      end do
      return
      end
