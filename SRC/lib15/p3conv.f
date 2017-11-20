      subroutine p3conv( itype, jrs, krs, jrf, krf )
c  Copyright (C) 2015, Jerry Sellwood
c
c performs the radial convolution part of the solution for the field
c   on the 3-D polar grid
c The Green function is presumed stored in Fourier transformed form
c   in a file (created by P3GRNM) and is read in to be used here.  If
c   the local array green is dimensioned to be large enough, it is
c   read in once at the first call and stored, otherwise it is read in
c   record-by-record as needed, at each call.
c The routine is called for each type of convolution, for radial forces,
c   azimuthal forces, vertical forces and potentials and the results
c   over-write the input data array
      use aarrays
      implicit none
c
c calling arguments
      integer itype, jrs, krs, jrf, krf
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
c local allocatable arrays
      real, allocatable, save :: gr(:)
c
      real, allocatable, save :: green(:,:,:,:)
c
      real, allocatable, save :: w(:)
c
c local variables
      integer g, i, ifail, im, ir, is, iuse, iz, j, jm, jz, k, km
      integer kz, l, m, mm, n, nit, nm, nrf, nrg, nrs, ntypes
      logical file, naskip, nvskip
      save file, iuse
      data iuse / 0 /
c
      ntypes = 4
      nrg = nr( jgrid )
      nit = ngt
c first call
      if( iuse .eq. 0 )then
        if( parallel )call crash( 'P3CONV',
     +       'Single processor version called for parallel operations' )
c check Green function file header
        call grnhed( .true., .false., ifail )
        if( ifail .ne. 0 )call crash( 'P3CONV',
     +                                     'Wrong Green function file' )
c determine whether green function is used from file and report
        file = .true.
c        file = .false.
        if( file )then
          if( master )then
            print *, 'Green fn used from file'
            write( no, * )'Green fn used from file'
          end if
        else
          if( master )then
            print *, 'Green function stored in memory'
            write( no, * )'Green function stored in memory'
          end if
c allocate permanent buffer for full Green fn
          l = nrg * nrg
          allocate ( green( l, ng, ngz + 1, ntypes ) )
          do j = 1, ntypes
            jm = 1
            km = ng
            if( j .eq. 2 )then
              jm = 2
              km = min( ng, mmax - 1 )
            end if
            jz = 1
            kz = ngz + 1
            if( j .eq. 3 )then
              jz = 2
              kz = ngz
            end if
            do iz = jz, kz
              do im = jm, km
                read( nit )( green( i, im, iz, j ), i = 1, l )
              end do
            end do
          end do
          close( nit )
        end if
c allocate local buffer for Green fn and work space
        l = nrg * nrg
        allocate ( gr( l ) )
        l = 2 * mesh( jgrid )
        allocate ( w( l ) )
        iuse = 1
      end if
c skip header if necessary
      if( file .and. ( itype .eq. 1 ) )then
        rewind nit
        read( nit )
      end if
c set sectoral harmonic range
      jm = 1
      km = ng
      if( itype .eq. 2 )then
        jm = 2
        km = min( ng, mmax - 1 )
      end if
      jz = 1
      kz = ngz + 1
      if( itype .eq. 3 )then
        jz = 2
        kz = ngz
      end if
c count active terms
      nm = 0
      do im = 1, ng
        if( .not. lg( im ) )then
          if( im .eq. 1 .or. im .eq. mmax )then
            nm = nm + 1
          else
            nm = nm + 2
          end if
        end if
      end do
c numbers of rings for sources and field
      nrs = krs - jrs + 1
      nrf = krf - jrf + 1
c clear results array
      k = 2 * nm * ngz * nrf
      do i = 1, k
        w( i ) = 0
      end do
c work over planes
      do iz = jz, kz
        nvskip = ( iz .gt. 1 ) .and. ( iz .le. ngz )
c work over sectoral harmonics
        mm = -1
        if( ( .not. lg( 1 ) ) .and. ( itype .eq. 2 ) )mm = 0
        do im = jm, km
c get Green function
          k = nrg * nrg
          if( file )then
            read( nit )( gr( g ), g = 1, k )
          else
            call blkcpy( green( 1, im, iz, itype ), gr, k )
          end if
c skip azimuthal wavenumber unless requested
          if( .not. lg( im ) )then
            naskip = ( im .gt. 1 ) .and. ( im .lt. mmax )
            mm = mm + 1
c m & n point to the cc and cs terms of the transformed mass array
            if( lg( 1 ) )then
              m = 1 + 2 * mm
            else
              m = max( 1, 2 * mm )
            end if
            m = ( m - 1 ) * 2 * ngz + max( 1, 2 * ( iz - 1 ) )
            n = m + 2 * ngz
c sum over radii for sources
            do is = jrs, krs
c j & k point to the cc and cs terms of the field - j + 1 & k + 1 to cs & ss
              if( lg( 1 ) )then
                j = 1 + 2 * mm
              else
                j = max( 1, 2 * mm )
              end if
              j = ( j - 1 ) * 2 * ngz + max( 1, 2 * ( iz - 1 ) )
              k = j + 2 * ngz
c g points to the Green function
              g = ( is - 1 ) * nrg + jrf
c work over radii for radial forces or potentials
              if( ( itype .eq. 1 ) .or. ( itype .eq. 4 ) )then
                do ir = jrf, krf
                  w( j ) = w( j ) + grdmss( m, 0 ) * gr( g )
                  if( naskip )w( k ) =
     +                     w( k ) + grdmss( n, 0 ) * gr( g )
                  if( nvskip )then
                    w( j + 1 ) =
     +                     w( j + 1 ) + grdmss( m + 1, 0 ) * gr( g )
                    if( naskip )w( k + 1 ) =
     +                     w( k + 1 ) + grdmss( n + 1, 0 ) * gr( g )
                  end if
                  j = j + 2 * ngz * nm
                  k = k + 2 * ngz * nm
                  g = g + 1
                end do
c work over radii for azimuthal forces (naskip always true)
              else if( itype .eq. 2 )then
                do ir = jrf, krf
                  w( j ) = w( j ) + grdmss( n, 0 ) * gr( g )
                  w( k ) = w( k ) - grdmss( m, 0 ) * gr( g )
                  if( nvskip )then
                    w( j + 1 ) =
     +                 w( j + 1 ) + grdmss( n + 1, 0 ) * gr( g )
                    w( k + 1 ) =
     +                 w( k + 1 ) - grdmss( m + 1, 0 ) * gr( g )
                  end if
                  j = j + 2 * ngz * nm
                  k = k + 2 * ngz * nm
                  g = g + 1
                end do
c work over radii for vertical forces (nvskip always true)
              else if( itype .eq. 3 )then
                do ir = jrf, krf
                  w( j ) = w( j ) + grdmss( m + 1, 0 ) * gr( g )
                  w( j + 1 ) =
     +                 w( j + 1 ) - grdmss( m, 0 ) * gr( g )
                  if( naskip )then
                    w( k ) =
     +                     w( k ) + grdmss( n + 1, 0 ) * gr( g )
                    w( k + 1 ) =
     +                 w( k + 1 ) - grdmss( n, 0 ) * gr( g )
                  end if
                  j = j + 2 * ngz * nm
                  k = k + 2 * ngz * nm
                  g = g + 1
                end do
              end if
              m = m + 2 * ngz * nm
              n = n + 2 * ngz * nm
c end loop over sources
            end do
c end loop over sectoral harmonics
          end if
        end do
c end loop over planes
      end do
c copy results back
      k = 2 * nm * ngz * nrf
      if( itype .lt. 4 )then
        do i = 1, k
          grdfld( i, itype ) = w( i )
        end do
      else
        do i = 1, k
          grdmss( i, 0 ) = w( i )
        end do
      end if
      return
      end
