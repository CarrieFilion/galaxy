      subroutine p3convp( itype, jrs, krs, jrf, krf )
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
      include 'mpif.h'
c
c local allocatable arrays
      real, allocatable, save :: gr(:)
c
      real, allocatable, save :: green( :, :, :, : )
c
      integer, allocatable, save :: req(:,:)
c
      integer, allocatable, save :: tag(:)
c
      real, allocatable, save :: w(:)
c
c local array
      integer stats( mpi_status_size )
c
c local variables
      integer g, i, ib, ifail, im, ir, is, iuse, iz, izz, j, jm
      integer jz, k, kbf, km, kz, l, m, mm, n, nb, nm
      integer nmpi, nrf, nrg, nrs, ntypes
      logical file, lskp( 4 ), naskip, nvskip
      save file, iuse
      data iuse / 0 /
c
!      if( master )print *, 'starting p3convp', jgrid, itype
      ntypes = 4
      nrg = nr( jgrid )
c first call
      l = nrg * nrg
      if( iuse .eq. 0 )then
        if( .not. parallel )call crash( 'P3CONVP',
     +       'Parallel version called for single processor operations' )
        if( lg( 1 ) )call crash( 'P3CONVP',
     +                    'Bug when m=0 term is skipped not yet fixed' )
c check Green function file header
        call grnhed( .true., .false., ifail )
        if( ifail .ne. 0 )call crash( 'P3CONV',
     +                                     'Wrong Green function file' )
c Green array should be held in memory for efficent parallel operation
c  but too much allocated memory causes problems
c        file = .false.
        file = .true.
        if( file )then
          if( master )then
            print *, 'Green function used from file'
            write( no, * ) 'Green function used from file'
          end if
        else
c read in Green function
          if( master )then
            print *, 'Green function stored in memory'
            write( no, * )'Green function stored in memory'
          end if
c allocate permanent buffer for full Green fn
          izz = ngz / numprocs + 1
          allocate ( green( l, ng, izz, ntypes ) )
          do j = 1, ntypes
            jm = 1
            km = ng
            if( j .eq. 2 )then
              jm = 2
              km = min( ng, mmax - 1 )
            end if
            jz = 1
            kz = ngz + 1
            if( fixrad .or. ( j .eq. 3 ) )then
              jz = 2
              kz = ngz
            end if
c read only planes needed by this node
            izz = 0
            do iz = jz, kz
              if( mod( iz - jz, numprocs ) .eq. myid )then
                izz = izz + 1
                do im = jm, km
                  read( ngt )( green( i, im, izz, j ), i = 1, l )
                end do
              else
c skip other planes
                do im = jm, km
                  read( ngt )
                end do
              end if
            end do
          end do
          close( ngt )
        end if
c allocate local buffer for Green fn and work space
        allocate ( gr( l ) )
        l = 2 * mesh( jgrid )
        allocate ( w( l ) )
        allocate ( req( numprocs, ngz + 1 ) )
        allocate ( tag( ngz + 1 ) )
        iuse = 1
        l = nrg * nrg
        if( master )print *, 'read in Green fn for grid', jgrid
      end if
c skip header if necessary
      if( file .and. ( itype .eq. 1 ) )then
        rewind ngt
        read( ngt )
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
      j = 0
      do im = jm, km
        if( .not. lg( im ) )j = j + 1
      end do
      nb = 4 * nrf * j
      k = nb * ( ngz + 1 )
c clear results array
      do i = 1, k
        w( i ) = 0
      end do
c maximum number of planes in mpi buffer
      nmpi = numprocs * 8
      nmpi = min( nmpi, 120 )
      do iz = 1, ngz + 1
        tag( iz ) = iz
      end do
c work over planes
      izz = 0
      do iz = jz, kz
c skip Green fn planes not needed on this node
        if( file )then
          if( mod( iz - jz, numprocs ) .ne. myid )then
            do im = jm, km
              read( ngt )
            end do
          end if
        end if
        if( mod( iz - jz, numprocs ) .eq. myid )then
          izz = izz + 1
          kbf = ( iz - 1 ) * nb
          nvskip = ( iz .gt. 1 ) .and. ( iz .le. ngz )
c work over sectoral harmonics
          mm = -1
          if( ( .not. lg( 1 ) ) .and. ( itype .eq. 2 ) )mm = 0
          do im = jm, km
c get Green function from file
            k = nrg * nrg
            if( file )read( ngt )( gr( g ), g = 1, k )
c skip azimuthal wavenumber unless requested
            if( .not. lg( im ) )then
              if( .not. file )call
     +                      blkcpy( green( 1, im, izz, itype ), gr, k )
              naskip = ( im .gt. 1 ) .and. ( im .lt. mmax )
              mm = mm + 1
c m & n point to the cc and cs terms of the transformed mass array
              if( lg( 1 ) )then
                m = 1 + 2 * mm
              else
                m = max( 1, 2 * mm )
              end if
              m = ( m - 1 ) * 2 * ngz + max( 1, 2 * ( iz - 1 ))
              n = m + 2 * ngz
c sum over radii for sources
              do is = jrs, krs
                ib = kbf
c g points to the Green function
                g = ( is - 1 ) * nrg + jrf
c work over radii for radial forces or potentials
                if( ( itype .eq. 1 ) .or. ( itype .eq. 4 ) )then
                  do ir = jrf, krf
                    w( ib + 1 ) = 
     +                     w( ib + 1 ) + grdmss( m, 0 ) * gr( g )
                    if( naskip )w( ib + 3 ) =
     +                     w( ib + 3 ) + grdmss( n, 0 ) * gr( g )
                    if( nvskip )then
                      w( ib + 2 ) =
     +                 w( ib + 2 ) + grdmss( m + 1, 0 ) * gr( g )
                      if( naskip )w( ib + 4 ) =
     +                 w( ib + 4 ) + grdmss( n + 1, 0 ) * gr( g )
                    end if
                    ib = ib + 4
                    g = g + 1
                  end do
c work over radii for azimuthal forces (naskip always true)
                else if( itype .eq. 2 )then
                  do ir = jrf, krf
                    w( ib + 1 ) =
     +                     w( ib + 1 ) + grdmss( n, 0 ) * gr( g )
                    w( ib + 3 ) =
     +                     w( ib + 3 ) - grdmss( m, 0 ) * gr( g )
                    if( nvskip )then
                      w( ib + 2 ) =
     +                 w( ib + 2 ) + grdmss( n + 1, 0 ) * gr( g )
                      w( ib + 4 ) =
     +                 w( ib + 4 ) - grdmss( m + 1, 0 ) * gr( g )
                    end if
                    ib = ib + 4
                    g = g + 1
                  end do
c work over radii for vertical forces (nvskip always true)
                else if( itype .eq. 3 )then
                  do ir = jrf, krf
                    w( ib + 1 ) =
     +                 w( ib + 1 ) + grdmss( m + 1, 0 ) * gr( g )
                    w( ib + 2 ) =
     +                     w( ib + 2 ) - grdmss( m, 0 ) * gr( g )
                    if( naskip )then
                      w( ib + 3 ) =
     +                 w( ib + 3 ) + grdmss( n + 1, 0 ) * gr( g )
                      w( ib + 4 ) =
     +                     w( ib + 4 ) - grdmss( n, 0 ) * gr( g )
                    end if
                    ib = ib + 4
                    g = g + 1
                  end do
                end if
c end loop over sources
                m = m + 2 * ngz * nm
                n = n + 2 * ngz * nm
              end do
c end loop over sectoral harmonics
              kbf = kbf + 4 * nrf
            end if
          end do
c initiate sends of this plane to all other nodes
          kbf = ( iz - 1 ) * nb + 1
          do i = 1, numprocs
            k = i - 1
            if( myid .ne. k )then
              call mpi_isend( w( kbf ), nb, mpi_real,
     +               k, tag( iz ), mpi_comm_world, req( i, iz ), ifail )
            end if
          end do
        else
c initiate receive of this plane by this node
          k = mod( iz - jz, numprocs )
          kbf = ( iz - 1 ) * nb + 1
          i = myid + 1
          call mpi_irecv( w( kbf ), nb, mpi_real, k, tag( iz ),
     +                             mpi_comm_world, req( i, iz ), ifail )
        end if
c wait for sends and receives from an earlier plane
        if( iz - nmpi .ge. jz )then
          if( mod( iz - nmpi - jz, numprocs ) .eq. myid )then
            do i = 1, numprocs
              j = iz - nmpi
              if( myid .ne. i - 1
     +                       )call mpi_wait( req( i, j ), stats, ifail )
            end do
          else
            i = myid + 1
            j = iz - nmpi
            call mpi_wait( req( i, j ), stats, ifail )
          end if
        end if
c end loop over planes
      end do
c clear the entire target msh array
      if( itype .lt. 4 )then
        do i = 1, 2 * nm * nrf * ngz
          grdfld( i, itype ) = 0
        end do
      else
        do i = 1, 2 * nm * nrf * ngz
          grdmss( i, 0 ) = 0
        end do
      end if
c wait for sends and receives from remaining planes
      izz = max( jz, kz - nmpi + 1 )
      do iz = izz, kz
        if( mod( iz - jz, numprocs ) .eq. myid )then
          do i = 1, numprocs
            k = i - 1
            if( myid .ne. k )call mpi_wait( req( i, iz ), stats, ifail )
          end do
        else
          i = myid + 1
          call mpi_wait( req( i, iz ), stats, ifail )
        end if
      end do
!      if( master )print *, 'done convolution'
c copy results back
      if( itype .lt. 4 )then
        do iz = jz, kz
          nvskip = ( iz .gt. 1 ) .and. ( iz .le. ngz )
          ib = ( iz - 1 ) * nb
          mm = -1
          if( ( .not. lg( 1 ) ) .and. ( itype .eq. 2 ) )mm = 0
          do im = jm, km
c skip azimuthal wavenumber unless requested
            if( .not. lg( im ) )then
              naskip = ( im .gt. 1 ) .and. ( im .lt. mmax )
              if( itype .eq. 1 )then
                lskp( 2 ) = nvskip
                lskp( 3 ) = naskip
                lskp( 4 ) = nvskip .and. naskip
              else if( itype .eq. 2 )then
                lskp( 2 ) = nvskip
                lskp( 3 ) = .true.
                lskp( 4 ) = nvskip
              else if( itype .eq. 3 )then
                lskp( 2 ) = .true.
                lskp( 3 ) = naskip
                lskp( 4 ) = naskip
              end if
              mm = mm + 1
c j & k point to the cc and cs terms of the field - j + 1 & k + 1 to cs & ss
              if( lg( 1 ) )then
                j = 1 + 2 * mm
              else
                j = max( 1, 2 * mm )
              end if
              j = ( j - 1 ) * 2 * ngz + max( 1, 2 * ( iz - 1 ))
              k = j + 2 * ngz
              do ir = jrf, krf
                grdfld( j, itype ) = w( ib + 1 )
                if( lskp( 2 ) )grdfld( j + 1, itype ) = w( ib + 2 )
                if( lskp( 3 ) )grdfld( k, itype ) = w( ib + 3 )
                if( lskp( 4 ) )grdfld( k + 1, itype ) = w( ib + 4 )
                ib = ib + 4
                j = j + 2 * ngz * nm
                k = k + 2 * ngz * nm
              end do
            end if
          end do
        end do
      else
        do iz = jz, kz
          nvskip = ( iz .gt. 1 ) .and. ( iz .le. ngz )
          ib = ( iz - 1 ) * nb
          mm = -1
          if( ( .not. lg( 1 ) ) .and. ( itype .eq. 2 ) )mm = 0
          do im = jm, km
c skip azimuthal wavenumber unless requested
            if( .not. lg( im ) )then
              naskip = ( im .gt. 1 ) .and. ( im .lt. mmax )
              lskp( 2 ) = nvskip
              lskp( 3 ) = naskip
              lskp( 4 ) = nvskip .and. naskip
              mm = mm + 1
c j & k point to the cc and cs terms of the field - j + 1 & k + 1 to cs & ss
              if( lg( 1 ) )then
                j = 1 + 2 * mm
              else
                j = max( 1, 2 * mm )
              end if
              j = ( ( jrf - 1 ) * nm + j - 1 ) * 2 * ngz +
     +                                          max( 1, 2 * ( iz - 1 ) )
              k = j + 2 * ngz
              do ir = jrf, krf
                grdmss( j, 0 ) = w( ib + 1 )
                if( lskp( 2 ) )grdmss( j + 1, 0 ) = w( ib + 2 )
                if( lskp( 3 ) )grdmss( k, 0 ) = w( ib + 3 )
                if( lskp( 4 ) )grdmss( k + 1, 0 ) = w( ib + 4 )
                ib = ib + 4
                j = j + 2 * ngz * nm
                k = k + 2 * ngz * nm
              end do
            end if
          end do
        end do
      end if
!      if( master )print *, 'finished p3convp', itype
      return
      end
