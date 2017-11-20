      subroutine arrest
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to re-center the model - ie, to move the center-of-mass to the
c   origin and optionally to reset the linear momentum to zero.
c This routine does not reassign zone numbers, though perhaps it should
c A call to this routine MUST be followed by calls to MASSET & FINDF before
c   the next time step begins
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      logical offgrd
      real tsfac
c
c local array
      real*8 tm( mcmp ), v( 3 ), x( 3 )
c
c local variables
      integer i, ipp, is, j, jst, mlists, n, offgrid
      logical off
      real*8 actmas, pmfac, tfac
c
c execute this only if requested
      if( netcom .or. netlmo )then
c assumes a single processor
        if( parallel )call crash( 'ARREST', 'mpi version needed' )
c initialise centre-of-mass arrays
        do i = 1, 6 * mcmp
          check( i ) = 0.
          if( i .le. ncmp )then
            nspop( i ) = 0
            tm( i ) = 0
          end if
        end do
c work through all lists of particles
        mlists = nlists
        if( .not. offanal )mlists = nlists - 1
        do ilist = 1, mlists
c find zone and grid code for this list
          call interpret
          n = 1
          if( parallel )n = myid + 1
          inext = islist( 1, ilist, n )
c work through groups of particles
          do while ( inext .ge. 0 )
            call gather( jst )
c work through group
            do is = 1, jst
              ipp = iflag( is )
c get time-step factor
              tfac = tsfac( iz( is ) )
c sum positions and momenta
              nspop( ipp ) = nspop( ipp ) + 1
              tm( ipp ) = tm( ipp ) + pwt( is )
              do i = 1, ndimen
                cm( i, ipp ) = cm( i, ipp ) + pwt( is ) * oldc( i, is )
                pp( i, ipp ) = pp( i, ipp ) +
     +                         pwt( is ) * oldc( i + ndimen, is ) / tfac
              end do
            end do
          end do
        end do
c find total mass of particles left on the grid
        pmfac = pmass / ( lscale**3 * ts**2 )
        actmas = 0
        n = 0
        do ipp = 1, ncmp
          n = n + nspop( ipp )
          actmas = actmas + tm( ipp )
        end do
        if( n .le. 0 )call crash( 'ARREST', 'Logical error' )
        actmas = actmas * pmfac
        if( pertbn )actmas = actmas + mptbr
c compute corrections
        do i = 1, ndimen
          x( i ) = 0
          v( i ) = 0
c no change to CoM position unless requested
          if( netcom )then
            do ipp = 1, ncmp
              x( i ) = x( i ) + pmfac * cm( i, ipp )
            end do
            if( pertbn )x( i ) = x( i ) + mptbr * xptrb( i ) * lscale
            x( i ) = x( i ) / actmas
          end if
c no change to velocity of CoM unless requested
          if( netlmo )then
            do ipp = 1, ncmp
              v( i ) = v( i ) + pmfac * pp( i, ipp )
            end do
            if( pertbn )v( i ) = v( i ) +
     +                                  mptbr * vptrb( i ) * lscale * ts
            v( i ) = v( i ) / actmas
          end if
        end do
c
c adjust positions and velocities of all particles
c
c initialize new linked list table
        n = myid + 1
        do i = 1, nlists
          islist( 2, i, n ) = -1
        end do
        offgrid = 0
c work through all lists of particles
        do ilist = 1, nlists
c find zone and grid code for this list
          call interpret
          inext = islist( 1, ilist, n )
c work through groups of particles
          do while ( inext .ge. 0 )
            call gather( jst )
c work through group
            do is = 1, jst
c get time-step factor
              tfac = tsfac( iz( is ) )
c adjust positions and momenta
              do i = 1, ndimen
                newc( i, is ) = oldc( i, is ) - x( i )
                newc( i + ndimen, is ) =
     +                            oldc( i + ndimen, is ) - v( i ) * tfac
              end do
c check for particles leaving the grid
              off = offgrd( is )
              if( ilist .lt. nlists )then
                offgrid = offgrid + 1
                if( off )call addoff( is )
              else
                label( is ) = igrd( iflag( is ) )
                if( .not. off )then
                  offgrid = offgrid - 1
                  noffm = noffm - 1
                  iz( is ) = nzones
                end if
              end if
            end do
c reset labels
            call relabl( jst )
c restore coordinates
            call scattr( jst )
          end do
        end do
        if( ( .not. offanal ) .and. ( offgrid .gt. 0 ) )then
          if( master )print *, offgrid,
     +                            ' particles pushed off grid in ARREST'
c          call crash( 'ARREST', 'Particles pushed off grid' )
        end if
c update permanent counters
        noff = noff + noffm
c        nshort = nshort + nshortm
c        do i = 1, 4
c          ncen( i ) = ncen( i ) + ncenm( i )
c        end do
c save new list origins
        n = myid + 1
        do i = 1, nlists
          islist( 1, i, n ) = islist( 2, i, n )
        end do
c flag old mass arrays as useless
        if( nzones .gt. 1 )then
          do i = 2, nzones
            lstep( 1, i ) = -1
            lstep( 2, i ) = -1
          end do
        end if
c update perturber coordinates
        if( pertbn )then
          do i = 1, 3
            xptrb( i ) = xptrb( i ) - x( i ) / lscale
            vptrb( i ) = vptrb( i ) - v( i ) * gvfac
          end do
        end if
c update xcen and pastp arrays
        do ipp = 1, ncmp
          n = igrd( ipp )
          do j = 1, npast
            ipast( j ) = ( 1 - j ) * icstp
            tfac = ts * real( ipast( j ) )
            do i = 1, 3
              xcen( i, n ) = xcen( i, n ) - x( i )
              pastp( i, j, n ) = xcen( i, n ) - tfac * v( i )
            end do
          end do
        end do
      end if
      return
      end
