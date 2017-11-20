      subroutine emngrp( jst )
c  Copyright (C) 2016, Jerry Sellwood
      use aarrays
      implicit none
c Adds contributions from current group of particles to the global integrals
c Called from ANLGRP with the time-centred coordinates already computed and
c   passed as the second argument
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
c externals
      real axipot, halpot, phcorr, satpot
c
c local allocatable array for time centered coordinates
      real, allocatable, save :: coords( :,: )
c
c local variables
      integer i, is, j, kgrd
      logical lfirst
      real E, r
      save lfirst
      data lfirst / .true. /
c
c allocate space
      if( lfirst )then
        allocate ( coords( 6, mbuff ) )
        lfirst = .false.
      end if
      if( lzones )then
c get time-centered coordindates
        call cencds( jst, coords )
      else
c coordinates not yet time centered or adjusted for zones (initialization call)
        do is = 1, jst
          do i = 1, 6
            coords( i, is ) = oldc( i, is )
          end do
        end do
      end if
c work through group
      do is = 1, jst
c ignore test particles
        if( .not. testp( iflag( is ) ) )then
          kgrd = label( is )
c compute particle's specific energy
          E = 0
          r = 0
          do i = 1, ndimen
            r = r + ( coords( i, is ) - xcen( i, kgrd ) )**2
            E = E + .5 * coords( i + ndimen, is )**2
          end do
          r = sqrt( r )
c exclude particles far from the current center
          if( r .lt. lscale )then
            E = E + gpot( is )
            if( pertbn )E = E + satpot( is )
c these options make no dynamical sense
            if( fixrad .or. suppl .or. rigidh )then
              r = 0
              do i = 1, ndimen
                r = r + coords( i, is )**2
              end do
              r = sqrt( r )
              if( fixrad )E = E + axipot( r )
              if( suppl  )E = E + phcorr( r )
              if( rigidh )E = E + halpot( r )
            end if
c is this a new low value?
            j = jebind( kgrd )
            if( E .lt. bindp( 5, j, kgrd ) )then
c replace the previously least bound
              bindp( 5, j, kgrd ) = E
              do i = 1, 3
                bindp( i, j, kgrd ) = coords( i, is )
              end do
              bindp( 4, j, kgrd ) = pwt( is )
c find and remember the new least bound
              E = bindp( 5, 1, kgrd )
              j = 1
              do i = 2, nebind
                if( bindp( 5, i, kgrd ) .gt. E )then
                  j = i
                  E = bindp( 5, i, kgrd )
                end if
              end do
              jebind( kgrd ) = j
            end if
          end if
        end if
      end do
      return
      end
