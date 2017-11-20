      subroutine phyprp( nact, lprop, prop, lglen, alspi, jlen, besc,
     +                   nLbin, nLa, lzman, zman, lzpr, zpbin,
     +                   lmoni, ehmon, nwring, wring )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Routine to reduce, and save into the .res file, data accumulated during an
c   analysis step.  Summary results are printed if requested.
c Called from MEASURE
c
c Data types processed are:
c   ANGM - the binned distribution of the z-component of angular momentum
c   INTG - mainly the global integrals of the system
c   LGSP - normalized results of the log spiral analysis
c   MONI - specific energies and z-angular momenta of a sample of particles
c   RNGS - rings of test particles
c   SPHB - normalized results of the spherical Bessel fn analysis
c   VFLD - binned means and dispersions of the disc velocity components
c   VFLH - binned means and dispersions of the halo velocity components
c   ZANL - normalized results of the sectorial bending analysis
c
c calling arguments
      integer jlen, lglen, lmoni, lprop, lzman, lzpr, nact, nLbin
      integer nwring
      integer nLa( nact, nLbin )
      real alspi( nact, 2, lglen ), besc( nact, 2, jlen )
      real ehmon( 6, lmoni ), prop( nact, lprop ), wring( nwring )
      real zman( nact, 2, lzman ), zpbin( nact, lzpr )
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
c local allocatable array
      real, allocatable :: w(:)
c
c local arrays
      character*12 a( 9 )
      character*4 bstr( 10 )
c
c local variables
      integer i, im, ip, ipp, j, k, m, n, nbin, nn
      real bhol, gvfac2, ketot, oke, petot, pmfac, pn, r1, r2, rk, t
      real te, tt, vol
      include 'inc/pi.f'
      equivalence ( rk, k )
c
      data a / 'mass in bin ', 'mean orbl vy', 'sigma v tngl',
     +         'mean rad vly', 'sigma v radl', '   mean z   ',
     +         '  sigma z   ', ' mean z vely', '  sigma v_z ' /
      data bstr / 'VFLD', 'ANGM', 'INTG', 'LGSP', 'SPHB', 'ZANL',
     +            'VFLH', 'MONI', 'ZPRF', 'RNGS' /
c
c gather global integrals data
      if( parallel )call crash( 'PHYPRP', 'mpi version needed' )
c integrals
      if( master )then
        gvfac2 = gvfac * gvfac
        pmfac = pmass / ( lscale**3 * ts**2 )
        do ipp = 1, ncmp
          n = nspop( ipp )
          if( n .gt. 0 )then
            do i = 1, 3
              cm( i, ipp ) = cm( i, ipp ) / ( lscale * popm( ipp ) )
              pp( i, ipp ) = pp( i, ipp ) * pmfac * gvfac
              for( i, ipp ) = for( i, ipp ) / ( lscale * ts**2 )
              ang( i, ipp ) = ang( i, ipp ) * pmfac / ( lscale**2 * ts )
            end do
            ke( ipp ) = .5 * ke( ipp ) * pmfac * gvfac2
            pe( ipp ) = pe( ipp ) * pmfac * gvfac2
            if( lprint )then
              write( no,
     +         '(/i10, '' particles analysed for icmp ='', i4 )' )n, ipp
              write( no, '( '' Total force components'', 3f10.5 )' )
     +             ( for( i, ipp ), i = 1, 3 )
              write( no, '( ''         Centre of mass'', 3f10.5 )' )
     +             ( cm( i, ipp ), i = 1, 3 )
              write( no, '( ''  Linear mom components'', 3f10.5 )' )
     +             ( pp( i, ipp ), i = 1, 3 )
              write( no, '( '' Angular mom components'', 3f10.5 )' )
     +             ( ang( i, ipp ), i = 1, 3 )
              write( no, '( '' Potential energy'', f10.4 )' )pe( ipp )
              write( no, '( ''   Kinetic energy'', f10.4 )' )ke( ipp )
            end if
          end if
        end do
        petot = selfe + pe( 2 ) + pe( 1 )
        ketot = ke( 1 ) + ke( 2 )
        te = petot + ketot
        claus = claus * pmfac * gvfac2
        if( lprint )then
          write( no, 200 )pe, selfe
 200      format( / '   Three pe terms are:', f10.4, ' Disc', f10.4,
     +            ' cross term and', f10.4, ' self energy of halo' )
          write( no, 201 )petot, ketot, claus
 201      format( 5x, 'Total pe', f10.4 / 5x, 'Total ke', f10.4 / 7x,
     +            'Virial', f10.4 )
          write( no, '( 2x, ''Total energy'', f10.4 )' )te
        end if
      end if
c normalize and save particle integrals
      if( moni )then
        if( master )then
          read( bstr( 8 ), '( a4 )' )bhol
          if( twod )then
            do i = 1, nmonit
              ehmon( 1, i ) = ehmon( 1, i ) * gvfac**2
              ehmon( 2, i ) = ehmon( 2, i ) * gvfac / lscale
              ehmon( 3, i ) = ehmon( 3, i ) / lscale
              ehmon( 4, i ) = ehmon( 4, i ) / lscale
              ehmon( 5, i ) = ehmon( 5, i ) * gvfac
              ehmon( 6, i ) = ehmon( 6, i ) * gvfac
            end do
            write( nphys )irun, bhol, istep, 6, nmonit
            write( nphys )( ( ehmon( j, i ), j = 1, 6 ), i = 1, nmonit )
          else
            do i = 1, nmonit
              ehmon( 1, i ) = ehmon( 1, i ) * gvfac**2
              ehmon( 2, i ) = ehmon( 2, i ) * gvfac / lscale
              ehmon( 3, i ) = ehmon( 3, i ) * gvfac / lscale
              ehmon( 4, i ) = ehmon( 4, i ) * gvfac / lscale
            end do
            write( nphys )irun, bhol, istep, 4, nmonit
            write( nphys )( ( ehmon( j, i ), j = 1, 4 ), i = 1, nmonit )
          end if
        end if
      end if
c save angular momentum data
      if( angm .and. master )then
        read( bstr( 2 ), '( a4 )' )bhol
        do ipp = 1, ncmp
          if( nspop( ipp ) .gt. 0 )then
            write( nphys )ipp, bhol, istep, nLbin, 1
            write( nphys )( nLa( ipp, i ), i = 1, nLbin )
            if( lprint )then
              write( no, * )
              write( no, * )'     Angular momentum distribution'
              write( no, * )
              write( no, '( 1x, 20i6 )' )( nLa( ipp, i ), i = 1, nLbin )
            end if
          end if
        end do
      end if
c normalize and save spherical Bessel coefficients
      if( sphb .and. uqmass )call crash( 'PHYPRP',
     +                        'Unequal particle masses not programmed' )
      if( master .and. sphb )then
        read( bstr( 5 ), '( a4 )' )bhol
        do ipp = 1, ncmp
          n = npbess( ipp )
          if( n .gt. 0 )then
            do i = 1, jlen
              besc( ipp, 1, i ) = besc( ipp, 1, i ) / real( n )
              besc( ipp, 2, i ) = besc( ipp, 2, i ) / real( n )
            end do
            write( nphys )ipp, bhol, istep, jlmax, jnmax
            write( nphys )( besc( ipp, 1, i ), i = 1, jlen ),
     +                    ( besc( ipp, 2, i ), i = 1, jlen )
          end if
        end do
      end if
c save rings data
      if( rngs .and. master )then
        read( bstr( 10 ), '( a4 )' )bhol
        write( nphys )irun, bhol, istep, nrings, npring
        write( nphys )( wring( i ), i = 1, nwring )
      end if
c
c normalize and save logarithmic spiral coefficients
      if( lgsp .and. master )then
        read( bstr( 4 ), '( a4 )' )bhol
        do ipp = 1, ncmp
          if( disc( ipp ) .and. ( nspop( ipp ) .gt. 0 ) )then
            k = 0
            do ip = 1, np
              do im = 1, nm
                k = k + 1
                alspi( ipp, 1, k ) = alspi( ipp, 1, k ) / popm( ipp )
                alspi( ipp, 2, k ) = alspi( ipp, 2, k ) / popm( ipp )
              end do
            end do
            write( nphys )ipp, bhol, istep, np, nm
            write( nphys )( alspi( ipp, 1, i ), i = 1, k ),
     +                    ( alspi( ipp, 2, i ), i = 1, k )
          end if
        end do
      end if
c normalize and save z mean analysis coefficients
      if( zanl .and. master )then
        read( bstr( 6 ), '( a4 )' )bhol
        do ipp = 1, ncmp
          if( disc( ipp ) .and. ( nspop( ipp ) .gt. 0 ) )then
            k = 0
            do ip = 1, nrz
              k = k + 1
              if( zman( ipp, 2, k ) .gt. 0. )then
                t = 1. / ( zman( ipp, 2, k ) * lscale )
                zman( ipp, 1, k ) = t * zman( ipp, 1, k )
                do im = 1, nmz
                  k = k + 1
                  zman( ipp, 1, k ) = t * zman( ipp, 1, k )
                  zman( ipp, 2, k ) = t * zman( ipp, 2, k )
                end do
              else
                k = k + nmz
              end if
            end do
            write( nphys )ipp, bhol, istep, nrz, nmz
            write( nphys )( zman( ipp, 1, i ), i = 1, k ),
     +                    ( zman( ipp, 2, i ), i = 1, k )
          end if
        end do
      end if
c normalize and save z-density profile of disk
      if( zprf .and. master )then
        read( bstr( 9 ), '( a4 )' )bhol
        do ipp = 1, ncmp
          if( disc( ipp ) .and. ( nspop( ipp ) .gt. 0 ) )then
            r2 = 0
            k = 0
            do i = 1, nbzr
              r1 = r2
              r2 = real( i ) / ( rbzf * lscale )
              vol = pi * ( r2**2 - r1**2 ) / ( zbzf * lscale )
              t = pmfac / vol
              do j = 1, nbzz
                k = k + 1
                zpbin( ipp, k ) = t * zpbin( ipp, k )
              end do
            end do
            write( nphys )ipp, bhol, istep, nbzr, nbzz
            write( nphys )( zpbin( ipp, i ), i = 1, lzpr )
          end if
        end do
      end if
c
c velocity field analyses
c
      if( ( vfld .or. vflh ) .and. master )then
        oke = 0.
        do ipp = 1, ncmp
c vertically averaged azimuthal analysis for discs
          if( vfld .and. disc( ipp ) .and. ( nspop( ipp ) .gt. 0 ) )then
            i = ndrad
            if( lprint )i = 3 * i
            allocate ( w( i ) )
c scan over measured properties
            do i = 1, nprop
              nn = 0
              do k = 1, ndrad
                w( k ) = 0.
              end do
c scan over angles
              do j = 1, ndang
                n = nn
c scan over radii
                do k = 1, ndrad
                  pn = prop( ipp, n + 1 )
                  if( pn .gt. 0. )then
                    m = n + i
                    tt = prop( ipp, m )
c compress data for print out
                    w( k ) = w( k ) + tt
c calculate averaged quantities
                    if( ( i .eq. 2 ) .or. ( i .eq. 4 )
     +                               .or. ( i .eq. 8 ) )then
                      if( i .eq. 2 )oke = oke + tt * tt / pn
                      tt = tt * gvfac / pn
                    else if( ( i .eq. 3 ) .or. ( i .eq. 5 )
     +                                    .or. ( i .eq. 9 ) )then
                      t = prop( ipp, m - 1 )
                      tt = tt * gvfac2
                      t = ( tt - pn * t * t ) / pn
                      t = max( t, 0. )
                      tt = sqrt( t )
                    else if( i .eq. 6 )then
                      tt = tt / ( lscale * pn )
                    else if( i .eq. 7 )then
                      t = prop( ipp, m - 1 )
                      tt = tt / ( lscale * lscale )
                      t = ( tt - pn * t * t ) / pn
                      t = max( t, 0. )
                      tt = sqrt( t )
                    end if
                    if( i .gt. 1 )prop( ipp, m ) = tt
                  end if
                  n = n + idskp
                end do
                nn = nn + nprop
              end do
c print out properties as a function of radius only
              if( lprint )then
                if( i .eq. 1 )then
                  do k = 1, ndrad
                    if( w( k ) .eq. 0. )go to 20
                  end do
                  k = ndrad + 1
   20             nbin = k - 1
                end if
                do k = 1, nbin
                  j = k + ndrad
                  n = j + ndrad
                  if( i .eq. 1 )then
                    w( n ) = w( k )
                  else if( ( i .eq. 2 ) .or. ( i .eq. 4 )
     +                                  .or. ( i .eq. 8 ) )then
                    w( k ) = w( k ) * gvfac / w( n )
                    w( j ) = w( k )
                  else if( ( i .eq. 3 ) .or. ( i .eq. 5 )
     +                                  .or. ( i .eq. 9 ) )then
                    t = w( j )
                    w( k ) = w( k ) * gvfac2
                    t = ( w( k ) - w( n ) * t * t ) / w( n )
                    t = max( t, 0. )
                    w( k ) = sqrt( t )
                  else if( i .eq. 6 )then
                    w( k ) = w( k ) / ( lscale * w( n ) )
                    w( j ) = w( k )
                  else if( i .eq. 7 )then
                    t = w( j )
                    w( k ) = w( k ) / ( lscale * lscale )
                    t = ( w( k ) - w( n ) * t * t ) / w( n )
                    t = max( t, 0. )
                    w( k ) = sqrt( t )
                  end if
                end do
                write( no, * )
                write( no, '( 50x, a12 /  )' )a( i )
                if( i .eq. 1 )then
                  write( no, '( 1x, 10f12.1 )' )( w( k ), k = 1, nbin )
                else
                  write( no, '( 1x, 10f12.4 )' )( w( k ), k = 1, nbin )
                end if
              end if
            end do
            n = idskp * ndrad
c save data to file
            read( bstr( 1 ), '( a4 )' )bhol
            write( nphys )ipp, bhol, istep, ndrad, ndang
            write( nphys )( prop( ipp, k ), k = 1, n )
c return local allocated space
            deallocate ( w )
c cylindrically averaged analysis for spheroidal population
          else if( vflh .and. ( .not. disc( ipp ) ) .and.
     +             ( nspop( ipp ) .gt. 0 ) )then
            i = nhrbin
            if( lprint )i = 3 * i
            allocate ( w( i ) )
c scan over measured properties
            do i = 1, nprop
              nn = 0
              do k = 1, nhrbin
                w( k ) = 0.
              end do
c scan over layers
              do j = 1, nhzbin
                n = nn
c scan over radii
                do k = 1, nhrbin
                  pn = prop( ipp, n + 1 )
                  if( pn .ne. 0. )then
                    m = n + i
                    tt = prop( ipp, m )
c compress data for print out
                    w( k ) = w( k ) + tt
c calculate averaged quantities
                    if( ( i .eq. 2 ) .or. ( i .eq. 4 )
     +                               .or. ( i .eq. 8 ) )then
                      if( i .eq. 2 )oke = oke + tt * tt / pn
                      tt = tt * gvfac / pn
                    else if( ( i .eq. 3 ) .or. ( i .eq. 5 )
     +                                    .or. ( i .eq. 9 ) )then
                      t = prop( ipp, m - 1 )
                      tt = tt * gvfac2
                      t = ( tt - pn * t * t ) / pn
                      t = max( t, 0. )
                      tt = sqrt( t )
                    else if( i .eq. 6 )then
                      tt = tt / ( lscale * pn )
                    else if( i .eq. 7 )then
                      t = prop( ipp, m - 1 )
                      tt = tt / ( lscale * lscale )
                      t = ( tt - pn * t * t ) / pn
                      t = max( t, 0. )
                      tt = sqrt( t )
                    end if
                    if( i .gt. 1 )prop( ipp, m ) = tt
                  end if
                  n = n + ihskp
                end do
                nn = nn + nprop
              end do
c print out properties as a function of radius only
              if( lprint )then
                if( i .eq. 1 )then
                  do k = 1, nhrbin
                    if( w( k ) .eq. 0. )go to 21
                  end do
                  k = nhrbin + 1
   21             nbin = k - 1
                end if
                do k = 1, nbin
                  j = k + nhrbin
                  n = j + nhrbin
                  if( i .eq. 1 )then
                    w( n ) = w( k )
                  else if( ( i .eq. 2 ) .or. ( i .eq. 4 )
     +                                  .or. ( i .eq. 8 ) )then
                    w( k ) = w( k ) * gvfac / w( n )
                    w( j ) = w( k )
                  else if( ( i .eq. 3 ) .or. ( i .eq. 5 )
     +                                  .or. ( i .eq. 9 ) )then
                    t = w( j )
                    w( k ) = w( k ) * gvfac2
                    t = ( w( k ) - w( n ) * t * t ) / w( n )
                    t = max( t, 0. )
                    w( k ) = sqrt( t )
                  else if( i .eq. 6 )then
                    w( k ) = w( k ) / ( lscale * w( n ) )
                    w( j ) = w( k )
                  else if( i .eq. 7 )then
                    t = w( j )
                    w( k ) = w( k ) / ( lscale * lscale )
                    t = ( w( k ) - w( n ) * t * t ) / w( n )
                    t = max( t, 0. )
                    w( k ) = sqrt( t )
                  end if
                end do
                write( no, * )
                write( no, '( 50x, a12 /  )' )a( i )
                if( i .eq. 1 )then
                  write( no, '( 1x, 10f12.1 )' )( w( k ), k = 1, nbin )
                else
                  write( no, '( 1x, 10f12.4 )' )( w( k ), k = 1, nbin )
                end if
              end if
            end do
            n = ihskp * nhrbin
c save data to file
            read( bstr( 7 ), '( a4 )' )bhol
            write( nphys )ipp, bhol, istep, nhrbin, nhzbin
            write( nphys )( prop( ipp, k ), k = 1, n )
c return local allocated space
            deallocate ( w )
          end if
        end do
c Ostriker's t parameter
        oke = .5 * oke * pmfac * gvfac2
        t = oke / abs( petot )
        if( lprint )write( no, '( '' Ostrikers t parameter'', f8.4 )' )t
      end if
c
      if( intg .and. master )then
        read( bstr( 3 ), '( a4 )' )bhol
        n = 18 * ncmp + 4
        write( nphys )irun, bhol, istep, n, ncmp
        j = 0
        allocate ( w( n ) )
        do ipp = 1, ncmp
          if( uqmass )then
            w( j + 1 ) = popm( ipp )
          else
            k = nspop( ipp )
            w( j + 1 ) = rk
          end if
          do i = 1, 3
            w( j + i + 1 ) = cm( i, ipp )
            w( j + i + 4 ) = pp( i, ipp )
            w( j + i + 7 ) = for( i, ipp )
            w( j + i + 10 ) = ang( i, ipp )
          end do
          w( j + 14 ) = pe( ipp )
          w( j + 15 ) = ke( ipp )
          j = j + 15
        end do
        write( nphys )te, sngl( claus ), t, noff, ( w( i ), i = 1, j ),
     +                ( ( angoff( i, ipp ), i = 1, 3 ), ipp = 1, ncmp )
c return local allocated space
        deallocate ( w )
      end if
      return
      end
