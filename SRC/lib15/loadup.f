      subroutine loadup
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to read particle positions and velocities from a file, convert
c   them into internal units and load them in memory as required by the
c   main galaxy code.
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      integer zoneno
      logical offgrd
c
c local arrays
      integer msp( 3 )
      real, allocatable :: inb(:)
c
c local variables
      integer discard, i, is, izon, j, jp, k, l, m, mcoor, nx
      logical lh, log
      real a, pm, t, wp
      real*8 xmass
c
      call gtreal( 'Enter time of .pcs file', t )
      istep = nint( t / ts )
c initialise / stfile /
      nx = 0
      is = 1
      j = myid + 1
      do i = 1, nlists
        islist( 2, i, j ) = -1
      end do
      discard = 0
c set flags that time step zones and time centering are not yet implemented
      lzones = .false.
      tcent = .false.
c open particle data file and read header
      call opnfil( ndistf, 'pcs', 'unformatted', 'old', 'seq', i )
      if( i .ne. 0 )call crash( 'LOADUP', 'Failed to open input file' )
      read( ndistf )msp, mcoor, l, t, pm, log
      if( master )print *, 'Checking header in LOADUP'
c check for expected information
      j = 0
      nbod = 0
      do i = 1, 3
        if( msp( i ) .ne. nsp( i ) )then
          if( master .and. ( j .eq. 0 ) )then
            j = j + 1
            print *, 'Expected', ( nsp( m ), m = 1, 3 )
            print *, '    Read', msp
          end if
c          call crash( 'LOADUP', 'Wrong number of particles' )
        end if
        nsp( i ) = msp( i )
        nbod = nbod + msp( i )
      end do
      if( numprocs .eq. 1 )then
        lpf = nbod * nwpp
      else
        lpf = nwpp * ( ( nbod - 1 ) / numprocs + 1 )
      end if
      if( lpf .gt. lptcls )call space( lptcls, lpf, 'ptcls', 'LOADUP' )
      wp = 1
      if( mod( mcoor, 2 ) .eq. 0 )then
        if( uqmass )call crash( 'LOADUP',
     +           'Individual particle massses, missing from .pcs file' )
        i = mcoor
      else
        if( .not. uqmass )call crash( 'LOADUP',
     + 'Individual particle massses should have been set in .dat file' )
        i = mcoor - 1
      end if
      if( i .ne. ncoor )call crash( 'LOADUP', 'Wrong dimension' )
      if( log )then
        if( .not. pertbn )call crash( 'LOADUP',
     +                        'Unexpected perturber data in .pcs file' )
        pertbn = .true.
        read( ndistf )xptrb, vptrb, mptbr, eptbr
      else
        if( pertbn )call crash( 'LOADUP',
     +                              'Perturber data missing from file' )
        pertbn = .false.
      end if
c set time step counter
      istep = nint( t / ts )
c set particle mass and infer active mass fractions
      pmass = pm * lscale**3 * ts**2
      do icmp = 1, ncmp
        if( ( .not. uqmass ) .and. ( nsp( icmp ) .gt. 0 ) )then
          if( cmpmas( icmp ) .gt. 0. )then
            fmass( icmp ) = pm * nsp( icmp ) / cmpmas( icmp )
          else
            fmass( icmp ) = 1
          end if
        end if
      end do
c allocate space for input buffer
      allocate( inb( l * mcoor ) )
      m = l
      jp = 0
c work over populations
      do icmp = 1, ncmp
        if( nsp( icmp ) .gt. 0 )then
        if( npr( icmp ) .le. 0 )npr( icmp ) = 1
        xmass = 0
        if( master )then
          if( disc( icmp ) )then
            if( cmpmas( icmp ) .gt. 0. )then
              write( no, '( // 19x, ''Setting up disc particles'' )' )
            else
              write( no, '( // 19x, ''Setting up ring particles'' )' )
            end if
          else
            write( no, '( // 19x, ''Setting up halo particles'' )' )
          end if
        end if
c allow for multiple grids
        call switch( igrd( icmp ) )
c determine whether this component contains the heavy anchor particle
        if( lheavy )then
          lh = .false.
          do i = 1, ngrid
            if( jebind( i ) .eq. icmp )lh = .true.
          end do
        end if
c process particles
        do k = 1, nsp( icmp )
c check space
          if( is .gt. mbuff )then
            if( nx .gt. lptcls )call crash( 'LOADUP',
     +                                        'Particle file overflow' )
            call relabl( mbuff )
            call scattr( mbuff )
            is = 1
          end if
c get next buffer of particles
          if( m .ge. l )then
            m = k - 1
            if( icmp .gt. 1 )then
              do i = 1, icmp - 1
                m = m + nsp( i )
              end do
            end if
            j = min( l, nbod - m ) * mcoor
            read( ndistf, iostat = izon )( inb( i ), i = 1, j )
            if( izon .ne. 0 )then
              print *, k, ' particles read for component', icmp
              print *, m, ' total so far', nbod, ' expected'
              call crash( 'LOADUP', 'Error reading .pcs file' )
            end if
            m = 0
            j = 0
          end if
c store scaled coordinates
          do i = 1, ndimen
            newc( i, is ) = inb( j + i ) * lscale
            newc( i + ndimen, is ) = inb( j + i + ndimen ) / gvfac
          end do
          if( uqmass )then
            wp = inb( j + mcoor )
            xmass = xmass + wp
          end if
          m = m + 1
          j = j + mcoor
c get zone
          if( nzones .gt. 1 .or. nguard .gt. 0 )then
            izon = zoneno( is, .true. )
          else
            izon = 1
          end if
c update counters
          if( izon .le. 0 )nshort = nshort + 1
          if( hole )then
            a = sqrt( newc( 1, is )**2 + newc( 2, is )**2 )
            if( a .lt. sngl( rhole ) )ncen( icmp ) = ncen( icmp ) + 1
          end if
c discard useless particles outside the old grid
          a = newc( 1, is )**2 + newc( 2, is )**2 + newc( 3, is )**2
          if( a .eq. 0. )then
            discard = discard + 1
            nsp( icmp ) = nsp( icmp ) - 1
            nbod = nbod - 1
c            lpf = lpf - nwpp
          else
c check for particle outside new grid
            if( offgrd( is ) )then
              izon = nlists
              noff = noff + 1
            end if
c assign labels
            iz( is ) = izon
            loc( is ) = nx
            iflag( is ) = icmp
            pwt( is ) = wp
            label( is ) = jgrid
c record location and node of heavy anchor particle - always the last
            if( centrd .and. lheavy .and. lh )then
              if( k .eq. nsp( icmp ) )then
                do i = 0, numprocs - 1
                  if( mod( jp, numprocs ) .eq. i )then
                    if( myid .ge. i )then
                      loch( jgrid ) = nx
                    else
                      loch( jgrid ) = nx - nwpp
                    end if
                    jebind( jgrid ) = i
                  end if
                end do
              end if
            end if
c increment counters when needed
            if( mod( jp, numprocs ) .eq. myid )then
              nx = nx + nwpp
              is = is + 1
            end if
            jp = jp + 1
          end if
        end do
c end main loop over populations
      end if
      if( uqmass .and. ( nsp( icmp ) .gt. 0 ) )fmass( icmp ) =
     +                                       xmass / dble( nsp( icmp ) )
      end do
      close( ndistf )
c clear any particles remaining in  / buffer /
      is = is - 1
      if( is .gt. 0 )then
        call relabl( is )
        if(
     +  nx .gt. lptcls )call crash( 'LOADUP', 'Particle file overflow' )
        call scattr( is )
      end if
c note number of particles outside grid, if any
      if( master )then
        if( discard .gt. 0 )then
          write( no, * )discard, ' particles discarded in LOADUP'
          print *, discard, ' particles discarded in LOADUP'
        end if
c note number of particles outside grid, if any
        if( noff .gt. 0 )then
          write( no, * )noff, ' particles outside the grid in LOADUP'
          print *, noff, ' particles outside the grid in LOADUP'
        end if
      end if
c update / stfile /
      j = myid + 1
      do i = 1, nlists
        islist( 1, i, j ) = islist( 2, i, j )
      end do
c scale all velocities for time step zones
      call adjvls( .true. )
c report new fmass
      if( master )then
        do icmp = 1, ncmp
          if( disc( icmp ) )then
            if( cmpmas( icmp ) .gt. 0. )print *,
     +                                'disk fmass set to', fmass( icmp )
          else
            if( cmpmas( icmp ) .gt. 0. )print *,
     +                                'halo fmass set to', fmass( icmp )
          end if
        end do
      end if
c adjust for new pmass
      call slfset
c initialize grid centering
      if( centrd )then
        if( istep .ne. 0
     +            )call crash( 'LOADUP', 'Centering not at the outset' )
        do jp = 1, ncmp
          l = igrd( jp )
c OK only if multiple pops on the same grid have the same CoM (usually true)
          do i = 1, 3
            xcen( i, l ) = lscale * comi( i, jp )
          end do
          do j = 1, npast
            ipast( j ) = ( 1 - j ) * icstp
            t = ipast( j )
            do i = 1, 3
              pastp( i, j, l ) = lscale * comi( i, jp ) +
     +                                     t * comi( i + 3, jp ) / gvfac
            end do
          end do
        end do
        call cenpth
      end if
      return
      end
