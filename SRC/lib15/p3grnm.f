      subroutine p3grnm( direct )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Routine to create Green function files for the 3-D polar grid
c   It creates tables of the functions for radial forces, azimuthal forces,
c   vertical forces and potentials, finding the Fourier tansform of each
c   and storing it in a file.
c It will also save the direct (untransformed) data in a separate file if
c   requested by the calling argument being set to .true.   The values in
c   the direct file will differ from those returned by the externals
c   (radfor, azifor & grnfun) if fewer than the maximum number of Fourier
c   components remains active.  This option is needed for the cross check
c   by direct convolution.
c If the files are already present, the routine simply reads the headers
c   and verifies that the values are consistent with the current run.
c
c This version uses single precision FFTPAK routines
c
c calling argument
      logical direct
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
      integer ltrig
      real, allocatable :: trig(:)
c
      integer lwork
      real, allocatable :: w(:)
c
c externals
      real azifor, grnfun, grofu, radfor, vrtfor
c
c local variables
      integer i, im, ir, is, it, itype, it1, it2, iz, iz1, iz2
      integer j, m, meshl, n, nl, nm, nrg, nz, ptclsl
      logical azimth, resyn, vert
      real rs, rf, strig, th, zf
c
      if( .not. p3d )call crash( 'P3GRNM', 'Wrong type of code' )
c open files
      if( ngt .lt. 0 )call new_unit( ngt )
      call opnfil( ngt, 'grt', 'unformatted', 'old', 'seq', i )
      if( i .eq. 0 )then
c verify existing transformed file
        call grnhed( .true., .false., i )
        if( i .ne. 0 )call crash( 'P3GRNM',
     +                             'Wrong transformed Green fun file' )
        if( direct )then
          if( ngd .lt. 0 )call new_unit( ngd )
          call opnfil( ngd, 'grd', 'unformatted', 'old', 'seq', i )
          if( i .eq. 1 )call crash( 'P3GRNM',
     +                               'Direct Green fn file not found' )
c verify existing direct file
          call grnhed( .true., .true., i )
          if( i .ne. 0 )call crash( 'P3GRNM',
     +                                  'Wrong direct Green fun file' )
        end if
        if( master )print *, ' Old files found and checked in P3GRNM'
c should not create the file(s) on the slave nodes as well
      else if( master )then
        nrg = nr( jgrid )
c preserve value lptcls
        ptclsl = lptcls
c allocate space
        meshl = nrg * ( na / 2 + 1 ) * ( ngz + 1 )
        allocate ( grdmss( meshl, 1 ) )
        lptcls = nrg * nrg * ( ngz + 1 ) * ng
        allocate ( ptcls( lptcls ) )
        lwork = 2 * nrg * ( ngz + 1 ) * ( na / 2 + 1 )
        allocate ( w( lwork ) )
        ltrig = 3 * max( na, ngz + 1 ) + 15
        allocate ( trig( ltrig ) )
c create new files
        print *, 'P3GRNM: Creating new Green fn files'
        call opnfil( ngt, 'grt', 'unformatted', 'new', 'seq', i )
        call grnhed( .false., .false., i )
        if( direct )then
          if( ngd .lt. 0 )call new_unit( ngd )
          call opnfil( ngd, 'grd', 'unformatted', 'new', 'seq', i )
          call grnhed( .false., .true., i )
        end if
c resynthesize unless all sectoral harmonics are included
        resyn = ng .lt. mmax
        do im = 1, ng
          resyn = resyn .or. lg( im )
        end do
c form 3 force components and then potentials
        itype = 0
        do while ( itype .lt. 4 )
          itype = itype + 1
c ignore azimuthal forces if ng = 1
          if( ( itype .eq. 2 ) .and. ( ng .eq. 1 ) )itype = itype + 1
          azimth = itype .eq. 2
          vert = itype .eq. 3
          it1 = 1
          it2 = mmax
          if( azimth )then
            it1 = 2
            it2 = mmax - 1
          end if
          iz1 = 1
          iz2 = ngz + 1
          if( vert )then
            iz1 = 2
            iz2 = ngz
          end if
          nz = iz2 - iz1 + 1
c scan over radii of sources
          nl = 1
          do is = 1, nrg
            m = 0
            rs = grofu( real( is - 1 ) )
c generate table of direct values
            do ir = 1, nrg
              rf = grofu( real( ir - 1 ) )
              do iz = iz1, iz2
                zf = dzg * real( iz - 1 )
                do it = it1, it2
                  th = alpha * real( it - 1 )
                  m = m + 1
                  if( itype .eq. 1 )grdmss( m, 1 ) =
     +                                          radfor( rs, rf, th, zf )
                  if( itype .eq. 2 )grdmss( m, 1 ) =
     +                                          azifor( rs, rf, th, zf )
                  if( itype .eq. 3 )grdmss( m, 1 ) =
     +                                          vrtfor( rs, rf, th, zf )
                  if( itype .eq. 4 )grdmss( m, 1 ) =
     +                                          grnfun( rs, rf, th, zf )
                end do
              end do
            end do
c save direct version if needed
            if( direct .and. ( .not. resyn ) )then
              write( ngd )( grdmss( i, 1 ), i = 1, m )
            end if
c initialize FFTPAK routines
            if( azimth )then
              nm = na / 2 - 1
              call rsinti( nm, trig )
            else
              nm = na / 2 + 1
              call rcosti( nm, trig )
            end if
            strig = 1. / real( na )
c Fourier analysis in angle - assemble data for FFT
            n = nrg * nz
            if( azimth )then
c find Fourier sine transform of azimuthal forces
              i = 1
              do ir = 1, n
                call rsint( nm, grdmss( i, 1 ), trig )
                i = i + nm
              end do
            else
c find Fourier cosine transform of radial or vertical forces or potentials
              i = 1
              do ir = 1, n
                call rcost( nm, grdmss( i, 1 ), trig )
                i = i + nm
              end do
            end if
c filter out unwanted harmonics
            m = 0
            do i = 1, n
              do it = it1, it2
                m = m + 1
                if( it .gt. ng )then
                  grdmss( m, 1 ) = 0
                else
                  if( lg( it ) )grdmss( m, 1 ) = 0
                end if
              end do
            end do
c re-synthesise direct version if needed
            if( direct .and. resyn )then
              do i = 1, m
                w( i ) = strig * grdmss( i, 1 )
              end do
              m = 1
              if( azimth )then
                do ir = 1, n
                  call rsint( nm, w( m ), trig )
                  m = m + nm
                end do
              else
                do ir = 1, n
                  call rcost( nm, w( m ), trig )
                  m = m + nm
                end do
              end if
c save new table of direct values
              m = nm * n
              write( ngd )( w( i ), i = 1, m )
            end if
c revise nm
            n = nm
            nm = ng - it1 + 1
            if( itype .eq. 2 )nm = min( ng, mmax - 1 ) - 1
c assemble array for Fourier transformation in z
            m = 0
            do ir = 1, nrg
              do im = 1, nm
                j = ( ir - 1 ) * n * nz + im
                do iz = iz1, iz2
                  m = m + 1
                  w( m ) = grdmss( j, 1 )
                  j = j + n
                end do
              end do
            end do
c initialize FFTs
            if( vert )then
              call rsinti( nz, trig )
            else
              call rcosti( nz, trig )
            end if
            strig = 1. / real( ngz )
c Fourier transformation in z
            n = nm * nrg
            m = 1
            if( vert )then
c find Fourier sine transform of vertical forces
              do i = 1, n
                call rsint( nz, w( m ), trig )
                m = m + nz
              end do
            else
c find Fourier cosine transform of radial or azimuthal forces or potentials
              do i = 1, n
                call rcost( nz, w( m ), trig )
                m = m + nz
              end do
            end if
c copy back
            m = 0
            n = nz * nm
            do ir = 1, nrg
              do iz = 1, nz
                j = ( ir - 1 ) * n + iz
                do im = 1, nm
                  m = m + 1
                  grdmss( m, 1 ) = w( j )
                  j = j + nz
                end do
              end do
            end do
c store transformed matrix
            call blkcpy( grdmss( 1, 1 ), ptcls( nl ), m )
            nl = nl + m
          end do
c reformat transformed version
          call p3gwrt( itype )
        end do
        deallocate ( ptcls )
        deallocate ( w )
        deallocate ( trig )
        deallocate ( grdmss )
c restore value of lptcls
        lptcls = ptclsl
c close and reopen files to ensure they are completely written
        close(ngt )
        call opnfil( ngt, 'grt', 'unformatted', 'old', 'seq', i )
        if( direct )then
          close( ngd )
          call opnfil( ngd, 'grd', 'unformatted', 'old', 'seq', i )
        end if
        write( no, '( '' Green function creation complete'' )' )
        print *, 'P3GRNM: Green fn files ready'
      end if
      return
      end
