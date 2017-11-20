      subroutine dumpp
c  Copyright (C) 2015, Jerry Sellwood
c
c Saves the information needed to resume the simulation from the current
c   moment.
c Many parameters are saved by a call to subroutine HEDREC, after which
c   this routine writes the:
c      starting points for the linked lists
c      table of supplementary forces and potentials (if needed)
c      most recent mass array (to enable a perfect restart)
c      old mass arrays (only if multiple time step zones are in use)
c   and finally
c      the current phase space positions of all the particles
c
c In parallel runs, the mass arrays and active particles on each
c   processor are saved in separate files
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
      integer i, ih, ix, izon, j, k, l, n
      parameter ( n = 5000 )
      logical savem
      real x
      equivalence ( ix, x )
c
c check that particles are ready to be dumped
      if( .not. tcent )call crash( 'DUMPP', 'tcent flag down' )
      if( .not. lzones )call crash( 'DUMPP', 'lzones flag down' )
c store header, pointer and halo values
      call opnfil( ndump, 'dmp', 'unformatted', 'unknown', 'seq', i )
      call hedrec( ndump, .false. )
      write( ndump )irun, istep, noff, ncen, nshort, nrsup, suppl, n,
     +                                   numprocs, angoff, amoff, rguard
      if( suppl )then
        if( .not. lsupst
     +      )call crash( 'DUMPP', 'Supplementary forces table not set' )
        write( ndump )hrfac,
     +                      ( ( htab( j, i ), j = 1, 3 ), i = 1, nrsup )
      end if
      if( pertbn )then
        if( exsphr .or. exgnrc )then
          write( ndump )pertbr
        else if( extbar )then
          write( ndump )omegab, bphase, accptb( 1 ),
     +             ( accpz( 1, 1, i ), accpz( 1, 2, i ), i = 1, nzones )
        else if( exspir )then
          write( ndump )omegsp, sphase, epspi
        else
          call crash( 'DUMPP', 'Unrecognized perturber' )
        end if
        if( nzones .gt. 1 )write( ndump )
     +                   ( lstep( 1, i ), lstep( 2, i ), i = 2, nzones )
      end if
      if( centrd )then
        if( lheavy )then
          write( ndump )icstp, kci, ipast, pastp,
     +                          ( jebind( i ), loch( i ), i = 1, ngrid )
        else
          write( ndump )icstp, kci, ipast, pastp
        end if
      end if
c count and output number of particles on this node
      k = 0
      do i = 1, nlists
        j = islist( 1, i, myid + 1 )
        do while ( j .ge. 0 )
          k = k + 1
          x = ptcls( j + nwpp )
          j = ix
        end do
      end do
      write( ndump )numprocs, myid, k, nwpp,
     +                       ( islist( 1, i, myid + 1 ), i = 1, nlists )
c save random generator information on each node
      call rvsave( ia, nia, xa, nxa, i, .false. )
      write( ndump )ia, xa
c save most recent mass array
      do ih = 1, ngrid
        call switch( ih )
        if( ncode .gt. 11 )call crash( 'DUMPP', 'Unrecognised method' )
        if( .not. ( noslfg .or. ldrct ) )then
          savem = master .or. p2d .or. p3d
          k = mesh( ih )
c save zone 1
          if( savem )then
            if( s3d )then
              write( ndump )( s3dfld( i, 1 ), i = 1, k )
            else if( sf2d )then
              write( ndump )( sfpfld( i, 1 ), i = 1, k ), newmaxr
            else if( sf3d )then
              write( ndump )( ( sfpfld( i, j ), i = 1, k ), j = 1, 2 )
            else if( scf )then
              write( ndump )( sfpfld( i, 1 ), i = 1, k )
            else if( lgrd )then
              write( ndump )( grdmss( i, 1 ), i = 1, k )
            else
              call crash( 'DUMPP', 'unrecognized method 1' )
            end if
c extra array for hybrid method
            if( hybrid .and. ( ih .gt. 1 ) )then
              if( s3d )then
                write( ndump )( s3dfld( i, 2 ), i = 1, k )
              else
                call crash( 'DUMPP', '2nd mass array not s3d' )
              end if
            end if
          end if
        end if
      end do
      call switch( 0 )
c save additional mass arrays (if needed)
      if( ( nzones .gt. 1 ) .and. ( .not. noslfg ) )then
        do ih = 1, ngrid
          call switch( ih )
          k = mesh( ih )
          savem = master .or. p2d .or. p3d
          if( savem )then
            do izon = 2, nzones
c old and new mass arrays
              if( s3d )then
                write( ndump )lstep( 1, izon ),
     +                      ( s3dmoz( i, izon, 1 ), i = 1, k )
                write( ndump )lstep( 2, izon ),
     +                      ( s3dmss( i, izon, 1 ), i = 1, k )
                if( hybrid .and. ( ih .eq. 2 ) )then
                  if( s3d )then
                    write( ndump )lstep( 1, izon ),
     +                      ( s3dmoz( i, izon, 2 ), i = 1, k )
                    write( ndump )lstep( 2, izon ),
     +                      ( s3dmss( i, izon, 2 ), i = 1, k )
                  else
                    call crash( 'DUMPP', '2nd grid not s3d' )
                  end if
                end if
              else if( lgrd )then
                write( ndump )lstep( 1, izon ),
     +                      ( grdmoz( i, izon ), i = 1, k )
                write( ndump )lstep( 2, izon ),
     +                      ( grdmss( i, izon ), i = 1, k )
              else if( lsfp )then
                write( ndump )lstep( 1, izon ),
     +                      ( sfpmoz( i, izon ), i = 1, k )
                write( ndump )lstep( 2, izon ),
     +                      ( sfpmss( i, izon ), i = 1, k )
              else
                call crash( 'DUMPP', 'Unrecognized method 2' )
              end if
            end do
          end if
        end do
        call switch( 0 )
      end if
c store particle data
      if( master )write( no,
     +          '( '' Dumping'', i10, '' particles in groups of'', i6 )'
     +                                                    )lpf / nwpp, n
c write out particle data
      l = n * nwpp
      k = 0
      do while ( k .lt. lpf )
        j = k + 1
        k = k + l
        k = min( k, lpf )
        write( ndump )( ptcls( i ), i = j, k )
      end do
c save offgrid particle index list if needed
      if( master )then
        if( offfor .and. ( noff .gt. 0 ) )write( ndump )
     +                                       ( loco( i ), i = 1, noff )
      end if
      close( unit = ndump )
      return
      end
