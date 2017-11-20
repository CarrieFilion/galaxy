      subroutine orbset
c  Copyright (C) 2015, Jerry Sellwood
c
c Main driving routine to set initial velocity components of only those
c   disc particles that depend on the actual gravitational field of the
c   active and rigid masses at the start.  These velocities may include
c   the orbital velocities for a pre-specified value of Toomre's Q
c   parameter, if these components were not previously set from a .dfn
c   file.  As the vertical balance of a thickened disc can be established
c   only after the gravitational field has been determined from the
c   grid, this routine still needs to be executed for thickened discs even
c   when the orbital motions were set from a DF.
c
c The magnitues of the required velocity dispersions, and asymmetric drift
c   if needed, are first tabulated by a call the QSETUP.  The individual
c   velocity components are assigned in SETORB by selecting values at
c   random from Gaussian distributions having these precalulated means and
c   dispersions in SETORB.   Note that a quiet start in position is
c   preserved by assigning equal velocities to all particles on a given ring.
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'inc/setup.f'
c
c local variables
      logical cold, lvset
      real mscale
c
      if( tcent )call crash( 'ORBSET', 'Coordinates are time centered' )
      if( lzones )call crash( 'ORBSET',
     +               'Velocities already adjusted for time step zones' )
c work over populations
      do icmp = 1, ncmp
        lvset = .false.
c skip if there is nothing to do
        if( nsp( icmp ) .gt. 0 ) then
          if( disc( icmp ) )then
c leave a 2-D disk with a DF alone
            if( master .and. dist( icmp ) )print *,
     +                     'Disk orbital velocities unchanged in orbset'
            if( .not. ( twod .and. dist( icmp ) ) )then
              cold = ( .not. dist( icmp ) ) .and.
     +               ( dfcns( 1, icmp ) .eq. 0. )
              lvset = .true.
            end if
          else
            cold = .false.
            lvset = cdft( icmp ) .eq. 'JEAN'
          end if
c tabulate required dispersions, both in-plane and vertically
          if( lvset )then
            call switch( igrd( icmp ) )
            lprint = .true.
            call qsetup
c set vels of disk particles
            bke = 0.
            ake = 0.
            call setorb( cold )
            mscale = pmass / ( lscale**3 * ts**2 )
            ake = .5 * real( npr( icmp ) ) * ake * mscale
            bke = .5 * real( npr( icmp ) ) * bke * mscale
            if( master .and. ( .not. cold ) )write( no,
     +               '( 3x, ''Kinetic energy before setorb'', f10.4,' //
     +                               ' '' and after'' f10.4 )' )bke, ake
            call switch( 0 )
          end if
        end if
      end do
      return
      end
