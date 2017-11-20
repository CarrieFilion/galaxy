      subroutine restrt
c  Copyright (C) 2015, Jerry Sellwood
c
c Restores the information from a previous dump to resume the simulation.
c Many parameters are read in by a call to subroutine HEDREC, after which
c   this routine reads back the:
c      starting points for the linked lists
c      table of supplementary forces and potentials (if needed)
c      most recent mass array (to enable a perfect restart)
c      old mass arrays (only if multiple time step zones are in use)
c   and finally
c      the current phase space positions of all the particles
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
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
c local arrays
      integer nia, nxa
      parameter ( nia = 9, nxa = 4 )
      integer ia( nia )
      real*8 xa( nxa )
c
c local variables
      integer i, ih, izon, j, k, l, n, nb
      real dt
c
      if( numprocs
     +         .gt. 1 )call crash( 'RESTRT', 'Parallel version needed' )
      call gtreal( 'Enter time of restrt files', dt )
      istep = nint( dt / ts )
c check header
      call opnfil( ndump, 'dmp', 'unformatted', 'old', 'seq', i )
      if( i .ne. 0 )call crash( 'RESTRT', 'No .dmp file found' )
      call hedrec( ndump, .true. )
c set self-force tables etc
      call slfset
c restore pointer and halo data
      read( ndump )irun, istep, noff, ncen, nshort, nrsup, suppl, n,
     +                                          j, angoff, amoff, rguard
      if( j .ne. numprocs )then
        if( master )print *,
     +        'numprocs is now', numprocs, ' but restart file is for', j
        call crash( 'RESTRT', 'numprocs mix up 1' )
      end if
c
      if( suppl )then
        read( ndump )hrfac,
     +                     ( ( htab( j, i ), j = 1, 3 ), i = 1, nrsup )
        lsupst = .true.
      end if
      if( pertbn )then
        if( exsphr .or. exgnrc )then
          read( ndump )pertbr
        else if( extbar )then
          read( ndump )omegab, bphase, accptb( 1 ),
     +             ( accpz( 1, 1, i ), accpz( 1, 2, i ), i = 1, nzones )
        else if( exspir )then
          read( ndump )omegsp, sphase, epspi
        else
          call crash( 'RESTRT', 'Unrecognized perturber' )
        end if
        if( nzones .gt. 1 )read( ndump )
     +                   ( lstep( 1, i ), lstep( 2, i ), i = 2, nzones )
      end if
c read moving grid centre information and initialize
      if( centrd )then
        if( lheavy )then
          read( ndump )icstp, kci, ipast, pastp,
     +                          ( jebind( i ), loch( i ), i = 1, ngrid )
        else
          read( ndump )icstp, kci, ipast, pastp
        end if
        call cenpth
c estimate position of grid center
        dt = real( istep - ipast( 1 ) ) / real( icstp )
        do j = 1, ngrid
          do i = 1, 3
            xcen( i, j ) = cenfit( 1, i, j ) * dt +
     +                     cenfit( 2, i, j ) * ( 1. - dt )
          end do
        end do
      end if
c read node specific information
      read( ndump )j, k, nb, l,
     +                       ( islist( 1, i, myid + 1 ), i = 1, nlists )
      if( j .ne. numprocs )then
        print *,
     +        'numprocs is now', numprocs, ' but restart file is for', j
        call crash( 'RESTRT', 'numprocs mix up 2' )
      end if
      if( k .ne. myid )then
        print *, 'myid is', myid, ' but restart file is for', k
        call crash( 'RESTRT', 'myid mix up' )
      end if
      if( l .ne. nwpp )then
        print *, 'nwpp is', nwpp, ' but restart file has', l
        call crash( 'RESTRT', 'nwpp mix up' )
      end if
      if( nb * l .gt. lptcls )then
        print *, 'number of particles in file',
     +                 nb, ' but there is space for only', lptcls / nwpp
        call crash( 'RESTRT', 'Particle number mix up' )
      end if
c restore random generator information
      read( ndump )ia, xa
      call rvsave( ia, nia, xa, nxa, i, .true. )
c restore most recent mass array
      do ih = 1, ngrid
        call switch( ih )
        if( ncode .gt. 11 )call crash( 'RESTRT', 'Unrecognised method' )
        if( .not. ( noslfg .or. ldrct ) )then
c read zone 1
          k = mesh( ih )
          if( s3d )then
            read( ndump )( s3dfld( i, 1 ), i = 1, k )
          else if( sf2d )then
            read( ndump )( sfpfld( i, 1 ), i = 1, k ), newmaxr
            maxr = min( newmaxr, maxmaxr )
          else if( sf3d )then
            read( ndump )( ( sfpfld( i, j ), i = 1, k ), j = 1, 2 )
          else if( scf )then
            read( ndump )( sfpfld( i, 1 ), i = 1, k )
          else if( lgrd )then
            read( ndump )( grdmss( i, 1 ), i = 1, k )
          else
            call crash( 'RESTRT', 'Unrecognized method 1' )
          end if
c extra array for hybrid method
          if( hybrid .and. ( ih .gt. 1 ) )then
            if( s3d )then
              read( ndump )( s3dfld( i, 2 ), i = 1, k )
            else
              call crash( 'RESTRT', '2nd mass array not s3d' )
            end if
          end if
        end if
c set jrad & krad
        if( p2d .or. p3d )then
          do i = 1, nzones
            jrad( i ) = 1
            krad( i ) = nr( ih )
          end do
        end if
      end do
      call switch( 0 )
c number of zones for next step
      do i = 1, nzones
        if( mod( istep, nstep( i ) ) .eq. 0 )mzone = i
      end do
c restore additional mass arrays (if needed)
      if( ( nzones .gt. 1 ) .and. ( .not. noslfg ) )then
        do ih = 1, ngrid
          call switch( ih )
          k = mesh( ih )
          do izon = 2, nzones
c old and new mass arrays
            if( s3d )then
              read( ndump )lstep( 1, izon ),
     +                      ( s3dmoz( i, izon, 1 ), i = 1, k )
              read( ndump )lstep( 2, izon ),
     +                      ( s3dmss( i, izon, 1 ), i = 1, k )
              if( hybrid .and. ( ih .eq. 2 ) )then
                if( s3d )then
                  read( ndump )lstep( 1, izon ),
     +                      ( s3dmoz( i, izon, 2 ), i = 1, k )
                  read( ndump )lstep( 2, izon ),
     +                      ( s3dmss( i, izon, 2 ), i = 1, k )
                else
                  call crash( 'RESTRT', '2nd grid not s3d' )
                end if
              end if
            else if( lgrd )then
              read( ndump )lstep( 1, izon ),
     +                      ( grdmoz( i, izon ), i = 1, k )
              read( ndump )lstep( 2, izon ),
     +                      ( grdmss( i, izon ), i = 1, k )
            else if( lsfp )then
              read( ndump )lstep( 1, izon ),
     +                      ( sfpmoz( i, izon ), i = 1, k )
              read( ndump )lstep( 2, izon ),
     +                      ( sfpmss( i, izon ), i = 1, k )
            else
              call crash( 'RESTRT', 'Unrecognized method 2' )
            end if
          end do
        end do
      end if
c
c recover particle data
      lpf = nb * nwpp
      if( lpf .gt. lptcls )call space( lptcls, lpf, 'ptcls', 'RESTRT' )
      l = n * nwpp
      k = 0
      do while ( k .lt. lpf )
        j = k + 1
        k = k + l
        k = min( k, lpf )
        read( ndump )( ptcls( i ), i = j, k )
      end do
c recover offgrid particle index list if present
      if( offfor .and. ( noff .gt. 0 ) )then
        if( master )read( ndump )( loco( i ), i = 1, noff )
      end if
      close( unit = ndump )
c set flags for state of particle coordinates
      tcent = .true.
      lzones = .true.
      return
      end