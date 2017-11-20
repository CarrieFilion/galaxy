      subroutine unload
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c Adjusts all particle positions and velocities to external units and
c    writes out the data in a semi-standard form
c
c If lurand is .T., the routine also shuffles the particles in each
c    population so that those in rings and/or selected from DF are
c    dispersed.  In this way, a uniform stride through the ptcls array
c    during execution of GALAXY selects a fair sample for particle plots
c    and monitoring.  This capability is useful only when creating the
c    initial pcs file, and should be turned off at other times.
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'mpif.h'
c
c local allocatable arrays
      real, allocatable :: outb(:)
      real, allocatable :: scr(:,:)
c
c external
      logical gtlogl
c
c local arrays
      integer tag( mnodes ), mb( mnodes ), stats( mpi_status_size )
c
c local variables
      integer i, ip, ipp, iproc, is, j, jst, k, kb, l, lb, lo, m, nb
      parameter ( l = 5000 )
      real pm, t
      equivalence ( pm, ip )
c
c time centre coordinates if needed
      if( tcent )call tmcent( .false. )
c undo velocity scaling for time step zones is needed
      if( lzones )call adjvls( .false. )
c rescale coordinates to external units
      iproc = 1
      if( parallel )iproc = myid + 1
      do ilist = 1, nlists
        islist( 2, ilist, iproc ) = -1
      end do
c work through all lists of particles
      do ilist = 1, nlists
c find zone and grid code for this list
        call interpret
        inext = islist( 1, ilist, iproc )
c work through groups of particles
        do while ( inext .ge. 0 )
          call gather( jst )
c erase useless coordinates of off-grid particles
          if( ( ilist .eq. nlists ) .and. ( .not. offanal ) )then
            do is = 1, jst
              do i = 1, ncoor
                newc( i, is ) = 0
              end do
            end do
          else
c scale to external units
            do is = 1, jst
              do i = 1, ndimen
                newc( i, is ) = oldc( i, is ) / lscale
                newc( i + ndimen, is ) = oldc( i + ndimen, is ) * gvfac
              end do
            end do
          end if
c save revised coordinates
          call scattr( jst )
        end do
      end do
c update list origin table
      do i = 1, nlists
        islist( 1, i, iproc ) = islist( 2, i, iproc )
      end do
c count number of particles on each node
      nb = 0
      do i = 1, nlists
        ip = islist( 1, i, iproc )
        do while ( ip .ge. 0 )
          nb = nb + 1
          pm = ptcls( ip + nwpp )
        end do
      end do
c gather this information to master node
      if( parallel )then
        call mpi_gather( nb, 1, mpi_integer, mb, 1, mpi_integer,
     +                   0, mpi_comm_world, i )
        if( master )then
c sanity check
          j = 0
          do i = 1, numprocs
            j = j + mb( i )
c ensure nb is the largest
            nb = max( nb, mb( i ) )
          end do
          if( j .ne. nbod )then
            print *, j, nbod
            call crash( 'UNLOAD', 'Particle count error' )
          end if
        end if
c set tags
        do i = 1, numprocs
          tag( i ) = i
        end do
      else
        mb( 1 ) = nb
      end if
c open a new output file
      if( master
     +     )call opnfil( ndistf, 'pcs', 'unformatted', 'new', 'seq', i )
      call mpi_bcast( i, 1, mpi_integer, 0,  mpi_comm_world, j )
      if( i .ne. 0 )then
        if( master )print *,'A file with this time stamp already exists'
        if( gtlogl( 'Do you want to overwrite it' ) ) then
          if( master )call opnfil(
     +                   ndistf, 'pcs', 'unformatted', 'old', 'seq', i )
          call mpi_bcast( i, 1, mpi_integer, 0,  mpi_comm_world, j )
        else
          call crash( 'UNLOAD', 'Declined to overwite pcs file' )
        end if
      end if
c write a short header
      t = ts * real( istep )
      k = ncoor
      if( uqmass )k = ncoor + 1
      pm = pmass / ( lscale**3 * ts**2 )
      if( master )then
        print *, 'Creating .pcs file at time', t
        write( ndistf )( nsp( i ), i = 1, 3 ), k, l, t, pm, pertbn
        if( pertbn )write( ndistf )xptrb, vptrb, mptbr, eptbr
      end if
c allocate the local buffers
      lb = ncoor
      if( uqmass )lb = ncoor + 1
      lb = l * lb
      allocate ( outb( lb ) )
      lo = l * nwpp
      allocate ( scr( lo, numprocs ) )
c set counters, m for overall total, is for # on each node, k for output buffer
      m = 0
      k = 0
      kb = 0
      is = 0
      ipp = 0
      do while ( m .lt. nbod )
c (re)fill the buffer
        if( is .ge. nb )then
          print *, is, nb, m
          call crash( 'UNLOAD', 'No more partcles to gather' )
        end if
        j = min( l, nb - is ) * nwpp
        call blkcpy( ptcls( is * nwpp + 1 ), scr( 1, myid + 1 ), j )
c copy these particles to the master node
        if( parallel )then
          i = 0
          do iproc = 1, numprocs - 1
            if( iproc .eq. myid )then
              call mpi_send( scr( 1, iproc + 1 ), j, mpi_real,
     +                          0, tag( iproc + 1 ), mpi_comm_world, i )
            else if( myid .eq. 0 )then
              call mpi_recv( scr( 1, iproc + 1 ), j, mpi_real,
     +               iproc, tag( iproc + 1 ), mpi_comm_world, stats, i )
            end if
          end do
        end if
c reset counter
        j = 0
c assemble output buffer and write out particle data on master node only
        if( master )then
          do while ( j .lt. lo )
            do iproc = 1, numprocs
              if( is .lt. mb( iproc ) )then
c output the buffer when full
                if( k .ge. lb )then
c shuffle if requested
                  if( lurand )then
                    i = ncoor
                    if( uqmass )i = i + 1
                    if( kb .eq. 0 )then
c all in one population
                      call shuffl( outb, l, i )
                    else
c avoid shuffling 2 populations together
                      ip = kb / i
                      call shuffl( outb, ip, i )
                      ip = l - ip
                      call shuffl( outb( kb + 1 ), ip, i )
                      kb = 0
                    end if
                  end if
                  write( ndistf )outb
                  k = 0
                end if
c counter for all particles
                m = m + 1
c check for new pop flag
                pm = scr( j + nwpp - 1, iproc )
                if( ip .ne. ipp )then
                  ipp = ip
                  print *, '1st particle in pop', ipp, ' is', m
                  kb = k
                end if
c copy the coordinates
                do i = 1, ncoor
                  outb( k + i ) = scr( j + i, iproc )
                end do
                k = k + ncoor
                if( uqmass )then
                  k = k + 1
                  outb( k ) = scr( j + ncoor + 1, iproc )
                end if
              end if
            end do
c counter for particles on each node and scratch buffer
            is = is + 1
            j = j + nwpp
          end do
        end if
        call mpi_bcast( is, 1, mpi_integer, 0, mpi_comm_world, i )
        call mpi_bcast( m, 1, mpi_integer, 0, mpi_comm_world, i )
      end do
c output any remaining particles
      if( master )then
        if( k .gt. 0 )then
c shuffle if requested
          if( lurand )then
            i = ncoor
            if( uqmass )i = i + 1
            if( kb .eq. 0 )then
c all in one population
              ip = k / i
              call shuffl( outb, ip, i )
            else
c avoid shuffling 2 populations together
              ip = kb / i
              call shuffl( outb, ip, i )
              ip = k / i - ip
              call shuffl( outb( kb + 1 ), ip, i )
            end if
          end if
          write( ndistf )( outb( i ), i = 1, k )
        end if
        close( ndistf )
      end if
      return
      end
