      subroutine interpret
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set the grid type and zone numbers of particles in the current
c   linked list.  The encryption rule has to be consistent with that used
c   by subroutine SCATTR:
c
c     ilist = izone + nguard + ( jlist - 1 ) * nzones
c
c with the exception
c     ilist = nlists   for particles off the grid
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/grids.f'
c
c check for nonsense
      if( ( ilist .lt. 1 ) .or. ( ilist .gt. nlists ) )then
        print *, 'ilist =', ilist
        call crash( 'INTERPRET', 'Impossible list number' )
      end if
c special case
      if( ilist .eq. nlists )then
        izone = nlists
        jlist = ngrid
c general case
      else
c determine the method type (jlist) for this list
        if( nplanes .gt. 1 )then
          jlist = ( ilist - 1 ) / nplanes + 1
          iplane = ilist - nplanes * ( jlist - 1 )
c force single plane for all except c3d
          if( igrid( jlist ) .ne. 2 )iplane = 1
        else
          jlist = ( ilist - 1 ) / ( nzones + nguard ) + 1
          iplane = 1
        end if
        if( ( jlist .le. 0 ) .or. ( jlist .gt. ngrid ) )then
          print *, ilist, jlist, iplane
          call crash( 'INTERPRET', 'Logical error' )
        end if
c determine the time step zone for this list
        izone = ilist - nguard - ( jlist - 1 ) * ( nzones + nguard )
        if( nplanes .gt. 1 )then
c force single time-step zone
          izone = 1
        else if( ( izone + nguard .le. 0 ) .or.
     +           ( izone .gt. nzones ) )then
          print *, ilist, jlist, izone - nguard, nzones
          call crash( 'INTERPRET', 'Impossible izone' )
        end if
      end if
      return
      end
