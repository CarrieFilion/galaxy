      subroutine bindel
c  Copyright (C) 2016, Jerry Sellwood
      use aarrays
c builds a sub-table of the gravitationally most bound particles of
c   of each mass component
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, j, jst, n
      logical local
c
      if( lbind )then
        local = lzones .and. ( .not. tcent )
        if( local )call adjvls( .false. )
c initialize search for most bound particles
        do j = 1, ngrid
          jebind( j ) = 1
          do i = 1, nebind
            bindp( 5, i, j ) = 100
          end do
        end do
c work through lists
        n = 1
        if( parallel )n = myid + 1
        do ilist = 1, nlists
c find zone and grid code for this list
          call interpret
          call switch( jlist )
c work through groups of particles
          inext = islist( 1, ilist, n )
          do while ( inext .ge. 0 )
            call gather( jst )
c get accelerations
            call getacc( jst )
c identify most bound particles
            call emngrp( jst )
          end do
        end do
c reset velocities for zones if they were undone above
        if( local )call adjvls( .true. )
      else
        call crash( 'BINDEL', 'Needless call' )
      end if
      return
      end
