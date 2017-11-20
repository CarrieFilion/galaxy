      subroutine step
c  Copyright (C) 2015, Jerry Sellwood
c
c Main driving routine for advancing all the particles forward by one step.
c   The acceleration components over the entire grid should have been
c   determined previously by a call to FINDF
c
c The particles are processed in groups, each group being worked through the
c   following sequence of calls:
c         GATHER:  collect the next group from the main particle storage area
c         GETACC:  look up acceleration components to be applied
c         STPGRP:  advance positions and velocities
c        (ANLGRP:  include contributions of these particles to measurements)
c         SCATTR:  replace new coordinates in the main particle storage area
c         MASSGN:  assign mass to the grid using new positions OR
c                  determine new maxr for SFP code using new positions
c
c The routine begins with calls to
c   MASSCL to zero out the mass array(s) in readiness for the new values
c   MEASURE to initialise the analysis software
c After all particles have been processed
c   STARCT checks the number of particles assigned to the grid
c   SCALED mulitplies the mass array(s) by the appropriate normalisation factor
c   (MEASURE is called again to complete the analysis and save the information)
c
c   MEASURE and ANLGRP are called only if this is an analysis step
c   REZGRP is needed only if multiple time steps are in use
c
c In the 3-D Cartesian code, the particles are sorted by planes of the grid to
c   enable all those in one plane to be processed before beginning the next.
c   This is so that acceleration components can be determined from the
c   potential more efficiently by a single call to DIFPOT for each plane.
c The particles in each zone are located through a linked list with the initial
c   location stored in the array ISLIST and the last particle has a zero
c   pointer.  A new linked list is created for the new positions, and the
c   initial locations are copied over the original set after all particles
c   have been processed.
      use aarrays
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
c allocatable arrays for analysis
      integer jlen, lglen, llz, lmoni, lprop, lzman, lzpr, nact, nbplot
      integer nwring
      integer nLbin
      integer, allocatable :: nLa(:,:)
      real, allocatable :: alspi(:,:,:)
      real, allocatable :: besc(:,:,:)
      real, allocatable :: dispy(:)
      real, allocatable :: ehmon(:,:)
      real, allocatable :: Lzval(:)
      real, allocatable :: prop(:,:)
      real, allocatable :: wring(:)
      real, allocatable :: zman(:,:,:)
      real, allocatable :: zpbin(:,:)
c
c local variables
      integer i, j, jst, n
      logical firsth
c
      if( parallel )call crash( 'STEP', 'Parallel version needed' )
      if( .not. tcent
     +            )call crash( 'STEP', 'Coordinates not time centered' )
      if( .not. lzones )call crash( 'STEP',
     +                   'Velocities not adjusted for time step zones' )
c initialize analysis if requested
      if( phys )then
c        nact = 0
c        do i = 1, ncmp
c          if( nsp( i ) .gt. 0 )nact = nact + 1
c        end do
        nact = ncmp
c select 10,000 particles per population
        plot = .false.
        if( ipict .gt. 0 )plot = mod( istep, ipict ) .eq. 0
        if( plot )then
          nbplot = 10000 * nact * ncoor
          allocate ( dispy( nbplot ) )
        end if
        if( vfld .or. vflh )then
          lprop = max( nhrbin * nhzbin, ndrad * ndang ) * nprop
          allocate ( prop( nact, lprop ) )
        end if
        if( angm )then
          nLbin = maxLz + 2
          allocate ( nLa( nact, nLbin ) )
        end if
        if( lgsp )then
          lglen = np * nm
          allocate ( alspi( nact, 2, lglen ) )
        end if
        if( lval )then
          llz = nbod
          allocate ( Lzval( nbod ) )
        end if
        if( moni )then
          lmoni = nmonit
          allocate ( ehmon( 6, lmoni ) )
        end if
        if( rngs )then
          nwring = ndimen * nrings * npring
          allocate ( wring( nwring ) )
        end if
        if( sphb )then
          jlen = jnmax * ( ( jlmax + 2 ) * ( jlmax + 1 ) ) / 2
          allocate ( besc( nact, 2, jlen ) )
        end if
        if( zanl )then
          lzman = nrz * ( nmz + 1 )
          allocate ( zman( nact, 2, lzman ) )
        end if
        if( zprf )then
          lzpr = nbzr * nbzz
          allocate ( zpbin( nact, lzpr ) )
        end if
        call measure( .false.,
     +                nact, lprop, prop, nbplot, dispy, lglen, alspi,
     +                jlen, besc, nLbin, nLa, lzman, zman, lzpr, zpbin,
     +                lmoni, ehmon, nwring, wring, llz, Lzval )
      end if
c number of zones for this step
      if( nzones .gt. 1 )then
        do i = 1, nzones
          if( mod( istep, nstep( i ) ) .eq. 0 )mzone = i
        end do
      else
        mzone = 1
      end if
c set mass arrays to zero
      call masscl( .true. )
c initialize temporary counters for this step
      noffm = 0
      nshortm = 0
      do i = 1, 4
        ncenm( i ) = 0
      end do
c initialize new or old lists - depending on which are to be stepped forward
      n = 1
      if( parallel )n = myid + 1
      do ilist = 1, nlists
        call interpret
        if( c3d .or. ( izone .le. mzone ) .or.
     +      ( offanal .and. ( mzone .eq. nzones ) .and.
     +                      ( ilist .eq. nlists ) ) )then
          islist( 2, ilist, n ) = -1
        else
          islist( 2, ilist, n ) = islist( 1, ilist, n )
          islist( 1, ilist, n ) = -1
        end if
      end do
      firsth = heavies
c work through all lists
      do ilist = 1, nlists
c find zone and grid code for this list
        call interpret
        call switch( jlist )
c update accelerations of heavies
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
c compute acc on perturber only if only one zone and ignore off-grid particles
        if( pertbn )sumptb = ( izone .le. 1 ) .and. ( nzones .eq. 1 )
        inext = islist( 1, ilist, n )
c work through groups of particles
        do while ( inext .ge. 0 )
          isubst = 1
          call gather( jst )
c get accelerations
          call getacc( jst )
c analysis
          if( phys )then
            call anlgrp( jst, nact, lprop, prop, nbplot, dispy,
     +                   lglen, alspi, jlen, besc, nLbin, nLa,
     +                   lzman, zman, lzpr, zpbin, lmoni, ehmon,
     +                   nwring, wring, llz, Lzval )
          end if
c step forward
          call stpgrp( jst )
c rezone
          call rezgrp( jst )
c particles moving with fractional time steps
          if( ( nguard .gt. 0 ) .and. ( izone .le. 1 ) )then
            call incare( jst )
          end if
c update labels etc
          call relabl( jst )
c store new coordinates
          call scattr( jst )
c compute accelerations on perturber if needed
          if( pertbn )then
            call pertbz( jst )
          end if
          if( .not. noslfg )then
c assign mass of particles - which may have switched grids
            if( heavies )then
              if( ilist .eq. nlists )then
c deal with any returning heavies first - duh! there cannot be any!
                jlist = igrd( 2 )
                call switch( jlist )
                mhyb = jgrid
                call massgn( jst )
                jlist = igrd( 1 )
              end if
              call switch( jlist )
              mhyb = jgrid
              call massgn( jst )
            else
              do j = jlist, ngrid
                call switch( j )
                mhyb = jgrid
                call massgn( jst )
              end do
              if( hybrid .and. ( jlist .eq. 1 ) )then
                call switch( 2 )
                mhyb = 1
                call massgn( jst )
                call switch( 1 )
              end if
            end if
          end if
        end do
      end do
c update list origin table
      do i = 1, nlists
        islist( 1, i, n ) = islist( 2, i, n )
      end do
      call switch( 0 )
c may need to work through particles a second time for some SFP methods
      if( sf2d .or. sf3d )call sfpsum
c update permanent counters
      noff = noff + noffm
      nshort = nshort + nshortm
      do i = 1, 4
        ncen( i ) = ncen( i ) + ncenm( i )
      end do
c complete analysis if requested
      if( phys )then
        call measure( .true.,
     +                nact, lprop, prop, nbplot, dispy, lglen, alspi,
     +                jlen, besc, nLbin, nLa, lzman, zman, lzpr, zpbin,
     +                lmoni, ehmon, nwring, wring, llz, Lzval )
        if( plot )deallocate ( dispy )
        if( vfld .or. vflh )deallocate ( prop )
        if( angm )deallocate ( nLa )
        if( lgsp )deallocate ( alspi )
        if( lval )deallocate ( Lzval )
        if( moni )deallocate ( ehmon )
        if( sphb )deallocate ( besc )
        if( rngs )deallocate ( wring )
        if( zanl )deallocate ( zman )
        if( zprf )deallocate ( zpbin )
      end if
c advance motion of perturber if present
      if( pertbn )then
        call ptbstp
      end if
c increment step number
      istep = istep + 1
c revise estimate of mass centroid when requested and regrid the mass
      if( centrd .and. ( mod( istep, icstp ) .eq. 0 ) )then
        call newxcn( .false. )
        call masset
      else
c combine masses from different zones
        call mascmb
c check particle count and rescale density distribution
        call ncheck
        call scaled
      end if
      return
      end
