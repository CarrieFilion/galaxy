      subroutine freqs
c  Copyright (C) 2015, Jerry Sellwood
c
c Tabulates, and stores if requested, the radial run of Omega, kappa & kappa_z
c   in the mid-plane at the current stage of evolution.  The values are
c   derived from the azimuthally averaged central or vertical attraction from
c   the combined active and rigid masses
      use aarrays
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
c local allocatable array
      real, allocatable :: w(:)
c
c externals
      real axipot, grofu, halpot, kappaz, meanrf, phcorr
c
c local variables
      character*4 bstr( 2 )
      integer i, ia, ir, j, k, l, m, nrad
      real ahol, fr, kappa, omega, rad, rmin, vcirc2, vgrad
c
      data bstr / 'FRQS', 'PNTL' /
c
      if( ( lprint .or. lsave ) .and. master )then
c make measurements on a sensible grid
        i = 1
        if( hybrid )i = ngrid
        call switch( i )
c save potentials if available and required
        if( potl .and. pntl )then
c allocate space
          allocate ( w( mesh( jgrid ) ) )
          read( bstr( 2 ), '( a4 )' )ahol
          l = irun
          if( hybrid )l = jgrid
          if( p2d .or. p3a )then
c save the entire grid
            m = 0
            if( p3a )then
              do ir = 1, nr( jgrid )
                rad = grofu( real( ir - 1 ) )
                do ia = 1, na
                  m = m + 1
                  w( m ) = grdfld( m, 3 )
c add external potentials
                  if( fixrad )then
                    w( m ) = w( m ) + axipot( rad )
                  else
                    if( suppl )w( m ) = w( m ) + phcorr( rad )
                    if( rigidh )w( m ) = w( m ) + halpot( rad )
                  end if
                  w( m ) = gvfac**2 * w( m )
                end do
              end do
            else
              m = mesh( jgrid )
              do i = 1, m
                w( i ) = gvfac**2 * grdfld( i, 3 )
              end do
            end if
            if( p2d )write( nphys )l, ahol, istep, nr( jgrid ), na
            if( p3a )write( nphys )l, ahol, istep, nr( jgrid ), ngz
            write( nphys )( w( i ), i = 1, m )
          else if( p3d )then
c extract values in the mid-plane
            j = na * ( ngz / 2 )
            m = 0
            do ir = 1, nr( jgrid )
              rad = grofu( real( ir - 1 ) )
              do ia = 1, na
                m = m + 1
                w( m ) = grdfld( j + ia, 4 )
c add external potentials
                if( fixrad )then
                  w( m ) = w( m ) + axipot( rad )
                else
                  if( suppl )w( m ) = w( m ) + phcorr( rad )
                  if( rigidh )w( m ) = w( m ) + halpot( rad )
                end if
                w( m ) = gvfac**2 * w( m )
              end do
              j = j + na * ngz
            end do
            write( nphys )l, ahol, istep, nr( jgrid ), na
            write( nphys )( w( i ), i = 1, m )
c return allocated space
            deallocate ( w )
          else
            call crash( 'FREQS', 'Potential array not programmed' )
          end if
        end if
c turn off accumulation of perturber accelerations
        if( pertbn )sumptb = .false.
        ilist = nlists - 1
c set array size
        rmin = 0
        if( p2d .or. p3d .or. sf2d .or. sf3d )rmin = rinh( jgrid )
        nrad = ( rgrid( jgrid ) - rmin ) / drfrqs
        allocate ( w( 4 * nrad ) )
c scan over radii measured in grid units
        vcirc2 = 1
        l = 0
        do while ( vcirc2 .gt. 0. )
          l = l + 1
          if( l .le. nrad )then
            rad = rmin + real( l ) * drfrqs
            fr = meanrf( rad )
c save average orbital velocity
            vcirc2 = rad * max( -fr, 0. )
            w( l ) = sqrt( vcirc2 ) * gvfac
          else
            vcirc2 = 0
          end if
        end do
        nrad = l - 1
        if( nrad .gt. 0 )then
c calculate Omega and kappa
          do i = 1, nrad
            omega = w( i ) * lscale / ( rmin + real( i ) * drfrqs )
            w( i + nrad ) = omega
            if( ( i .eq. 1 ) .or. ( i .eq. nrad ) ) then
              w( i + 2 * nrad ) = 0.
            else
              vgrad = ( w( i + 1 ) - w( i - 1 ) ) * lscale /
     +                                                   ( 2. * drfrqs )
              kappa = 2. * omega * ( omega + vgrad )
              kappa = max( kappa, 0. )
              kappa = sqrt( kappa )
              w( i + 2 * nrad ) = kappa
            end if
            if( threed )then
              rad = ( rmin + real( i ) * drfrqs ) / lscale
              w( i + 3 * nrad ) = kappaz( rad )
            end if
          end do
c write out
          m = 3
          if( threed )m = 4
          if( lprint )then
            k = 0
            do l = 1, m
              if( l .eq. 1 )write( no, * )'0Local circular velocities'
              if( l .eq. 2 )write( no, * )'0Mean angular frequencies'
              if( l .eq. 2 )write( no, '( / 50x, ''omega'' / )' )
              if( l .eq. 3 )write( no, '( / 50x, ''kappa'' / )' )
              if( l .eq. 4 )write( no, '( / 50x, ''kappa_z'' / )' )
              j = k + 1
              k = k + nrad
              write( no, '( 10f12.5 )' )( w( i ), i = j, k )
            end do
          end if
c save frequencies if requested - ignore circular vels
          if( frqs )then
            read( bstr( 1 ), '( a4 )' )ahol
            write( nphys )irun, ahol, istep, nrad, m - 1
            j = nrad + 1
            k = m * nrad
            write( nphys )( w( i ), i = j, k )
          end if
        else
          print *,
     +           'warning: freqs requested, but none are being recorded'
        end if
        deallocate ( w )
c restore accumulation of perturber accelerations
        if( pertbn )sumptb = .true.
      end if
      return
      end
