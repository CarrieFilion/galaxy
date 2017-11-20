      subroutine kecmps
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to plot kinetic energy components from the runXXX.res file
c    the data for disc particles have been binned differently from those
c    for halo particles
c There are 5 types of data for 2-D discs and 9 for 3-D discs or halos
c
c Graphics routines are from JSLIB
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c external
      real roundup
c
c local arrays
      real, allocatable :: ke( : )
      real, allocatable :: ket( : )
      real, allocatable :: t( : )
c
c local variables
      character type*4
      integer i, ifail, iopt, j, l, ll, m, mt, nc, ndsc, nhal
      logical disk, halo
      real fnp, vm, vm2, vz, w, xmax, xmin, y, ymax, ymin
      equivalence ( m, w )
c
c determine which population to plot
      if( twod )then
        disk = .true.
        halo = .false.
        type = 'VFLD'
      else
        disk = nsp( 1 ) .gt. 0
        halo = nsp( 2 ) .gt. 0
        if( disk .and. halo )then
          disk = .false.
          halo = .false.
          do while ( .not. ( disk .or. halo ) )
            call gtchar( 'Look at disc or halo velocity field?', type )
            disk = type .eq. 'disc'
            halo = type .eq. 'halo'
          end do
        else
          if( .not. ( disk .or. halo ) )call crash( 'KECMPS',
     +                             'No particles in either population' )
        end if
        if( disk )type = 'VFLD'
        if( halo )type = 'VFLH'
      end if
c
c set scaling constants
      if( ( .not. disc( 1 ) ) .or. disc( 2 ) )call crash( 'KECMPS',
     +                                     'Pops not of expected form' )
      if( halo )then
        fnp = cmpmas( 2 ) * fmass( 2 ) / real( nsp( 2 ) )
      else
        fnp = cmpmas( 1 ) * fmass( 1 ) / real( nsp( 1 ) )
      end if
c allocate space
      mt = 2501
      allocate ( ke( mt ) )
      allocate ( ket( mt ) )
      allocate ( t( mt ) )
c select option
    1 call gtintg( 'Enter option, -1 to stop or 0 for help', iopt )
      if( iopt .lt. 0 )return
c primitive help menu
      if( iopt .eq. 0 )then
        print *, 'Options are:'
        print *, ' 1 - KE_a'
        print *, ' 2 - KE_r'
        if( threed )then
          print *, ' 3 - KE_z'
        end if
        go to 1
      end if
      if( twod .and. ( iopt .gt. 2 ) )iopt = iopt + 1
      if( iopt .gt. 3 )go to 1
c convert iopt to a pointer
      iopt = 1 + 2 * iopt
      if( iopt .eq. 7 )iopt = 9
c
c restart file
      call rewnd
      nc = 0
      ymax = 0
      call nextrec( type, ifail )
c start of main loop
      do while ( ifail .eq. 0 )
c work over radii
        vm = 0
        vm2 = 0
c sum over radii & angles or planes
        do i = 1, mr
          l = ( i - 1 ) * nprop * ma - nprop
          do j = 1, ma
            l = l + nprop
            ll = l + iopt
            vm = vm + wres( l + 1 ) * wres( ll - 1 )
            vm2 = vm2 + wres( l + 1 ) * ( wres( ll ) * wres( ll ) +
     +                                 wres( ll - 1 ) * wres( ll - 1 ) )
          end do
        end do
        call nextrec( 'INTG', ifail )
c take out bulk motion of galaxy for vertical component only
        if( iopt .eq. 9 )then
          w = wres( 32 )
          ndsc = m
          w = wres( 33 )
          nhal = m
          vz = wres( 9 ) + wres( 12 )
          i = ndsc
          if( halo )i = nhal
          vm2 = vm2 - 2. * vz * vm + real( i ) * vz**2
        end if
c save result
        nc = nc + 1
        if( nc .gt. mt )call crash( 'KECMPS', 'Work arrays too small' )
        t( nc ) = time
        ke( nc ) = .5 * fnp * vm2
c        i = 27
c        if( halo )i = 28
c        ket( nc ) = wres( i )
        ymax = max( ymax, ke( nc ), ket( nc ) )
c get next record
        call nextrec( type, ifail )
      end do
c set frame
      call jspage
      call jssize( .2, .95, .1, .9 )
      xmin = 0
      xmax = roundup( .99 * t( nc ) )
      ymin = 0
      ymax = roundup( ymax )
      call jscale( xmin, xmax, ymin, ymax )
      call jsaxis( 'x', 'time', 1 )
      if( iopt .eq. 3 )call jsaxis( 'y', 'KE_{\\phi}', 1 )
      if( iopt .eq. 5 )call jsaxis( 'y', 'KE_{r}', 1 )
      if( iopt .eq. 9 )call jsaxis( 'y', 'KE_{z}', 1 )
      call jsbldt( 'Run no' )
      call jsbldi( irun, 5 )
      if( disk )call jsbldt( 'disc' )
      if( halo )call jsbldt( 'halo' )
      y = 1.05 * ymax - .05 * ymin
      call jswrit( .95 * xmin + .05 * xmax, y )
c plot data
      call jsjoin( t, ke, nc )
c      call jsdash( 0, 2, 0, 2 )
c      call jsjoin( t, ket, nc )
c      call jsdash( 0, 0, 0, 0 )
      go to 1
      end
