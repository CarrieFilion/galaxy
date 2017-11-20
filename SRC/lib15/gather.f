      subroutine gather( jst )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Routine to gather the next group of particles belonging to the current
c   list into the arrays in / buffer /
c The particles have the appropriate number of phase space coordinates,
c   ncoor (= 4 or 6) stored in floating point plus an additional two
c   integers.
c The first is used as a flag, which contains the population to which the
c   particle belongs.  For 2-D models only so far, it might also contain
c   an encoded length of time step for particles on extra short steps
c   (zone=0) and/or a tag if this particle is selected for monitoring of
c   its integrals etc.
c The second is the link to the next particle in the list (or a negative
c   value if the list ends with this particle)
c
c calling arguments - input value of jst is ignored, returned value is
c   the actual number of particles gathered.  The maximum is set by the
c   parameter mbuff which is used to dimension the workspace arrays
c
c calling argument
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
c local variables
      integer i, ia
      real a
      equivalence ( ia, a )
c
      jst = 0
c gather particles from linked list
      do while ( jst .lt. mbuff )
        jst = jst + 1
c set unstored info
        ncl( jst ) = 0
        loc( jst ) = inext
        iz( jst ) = izone
        label( jst ) = jlist
        ipln( jst ) = iplane
c pick up coordinates
        do i = 1, ncoor
          oldc( i, jst ) = ptcls( inext + i )
        end do
c particle weight, if not uniform
        if( uqmass )pwt( jst ) = ptcls( inext + ncoor + 1 )
c Barnes-Hut tree loc
        if( bht )then
          a = ptcls( inext + nwpp - 2 )
          ltree( jst ) = ia
        end if
c pick up flag
        a = ptcls( inext + nwpp - 1 )
        iflag( jst ) = ia
c pick up next pointer
        a = ptcls( inext + nwpp )
        inext = ia
c check for end of list
        if( inext .lt. 0 )return
      end do
      return
      end
