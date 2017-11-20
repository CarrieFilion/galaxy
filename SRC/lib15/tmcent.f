      subroutine tmcent( start )
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to initialize or undo time centering of the coordinates
c If the logical variable start is .TRUE. then time centering is initialized
c   else it is undone.
c If LFRV is selected, the velocities only are shifted half a step backwards
c   (forwards) in order to set up (undo) time centering.  The positions are
c   unaffected.
c If APLF is selected, the positions only are shifted half a step forwards
c   (backwards) in order to set up (undo) time centering.  The velocities
c   are unaffected.
c For LFRV, this routine should be preceded by a call to findf in order that
c   the correct acceleration components can be determined in the usual way.
c The routine process all the particles.
      use aarrays
      implicit none
c
c calling argument
      logical start
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c local variables
      integer i, j, jst, n
      logical firsth
      real as
c
c check that this is a sensible call
      if( .not.
     +     lzones )call crash( 'TMCENT', 'Vels not adjusted for zones' )
      if( start .eqv. tcent )then
        if( start )then
          if( master )print *, 'tcent flag already up'
        else
          if( master )print *, 'tcent flag already down'
        end if
        call crash( 'TMCENT', 'Logic error in call' )
      end if
      if( master )then
        if( start )then
          print *, 'Time-centering velocities'
          write( no, * )'Time-centering velocities'
        else
          print *, 'Un-centering velocities'
          write( no, * )'Un-centering velocities'
        end if
      end if
c initialize new linked list table
      n = 1
      if( parallel )n = myid + 1
      do i = 1, nlists
        islist( 2, i, n ) = -1
      end do
      firsth = heavies
c work through all lists of particles
      do ilist = 1, nlists
c find zone and grid code for this list
        call interpret
        call switch( jlist )
c update accelerations of heavies
        if( ( lfrv ) )then
          if( firsth .and. ( jlist .gt. 1 ) )then
            do i = 1, ndrct
              do j = 1, 4
                drpt( j + 7, i ) = drpt( j + 7, i ) + drpt( j + 11, i )
              end do
            end do
            firsth = .false.
            stphev = .true.
          else
            if( heavies )stphev = .false.
          end if
        end if
        inext = islist( 1, ilist, n )
c work through groups of particles
        do while ( inext .ge. 0 )
          isubst = 1
          call gather( jst )
c get acceleration components only if needed
          if( lfrv )then
            call getacc( jst )
          end if
c time centre velocities
          call tmcgrp( jst, start )
c resave coordinates
          if( .not. lfrv )call relabl( jst )
          call scattr( jst )
        end do
      end do
c update permanent counters
      noff = noff + noffm
c      nshort = nshort + nshortm
c      do i = 1, 4
c        ncen( i ) = ncen( i ) + ncenm( i )
c      end do
c save new list origins
      do i = 1, nlists
        islist( 1, i, n ) = islist( 2, i, n )
      end do
c flag old mass arrays as useless - this if should never be true!
      if( ( .not. lfrv ) .and. ( nzones .gt. 1 ) )then
        do i = 2, nzones
          lstep( 1, i ) = -1
          lstep( 2, i ) = -1
        end do
      end if
c time centre motion of perturber - if present
      if( pertbn )then
        as = .5 * ts
        if( start )as = -as
        if( exsphr .or. exgnrc )then
          if( lfrv )then
            do i = 1, 3
              vptrb( i ) = vptrb( i ) + as * accptb( i )
            end do
          else
            do i = 1, 3
              xptrb( i ) = xptrb( i ) - as * vptrb( i )
            end do
          end if
        else if( extbar )then
          if( lfrv )then
            omegab = omegab + as * accptb( 1 )
          else
            bphase = bphase - as * omegab
          end if
        else if( exspir )then
          if( .not. lfrv )sphase = sphase - as * omegsp
        else
          call crash( 'TMCENT', 'Unrecognized perturber' )
        end if
      end if
c reset flag to indicate routine was executed
      tcent = start
      return
      end
