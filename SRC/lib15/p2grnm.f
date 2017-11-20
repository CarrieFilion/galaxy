      subroutine p2grnm( direct )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Routine to create Green function files for the 2-D polar grid
c   It creates tables of the functions for radial forces, azimuthal forces
c   and potentials, finding the Fourier transform of each and storing it
c   in a file.
c It will also save the direct (untransformed) data in a separate file if
c   requested by the calling argument being set to .true.   The values in
c   the direct file will differ from those returned by the the externals
c   (radfor, azifor & grnfun) if fewer than the maximum number of Fourier
c   components remains active.  This option is needed for the cross check
c   by direct convolution, and could also be used for potential energy
c   measurements, though it is not at present.
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
      real azifor, grnfun, grofu, radfor
c
c local variables
      integer i, ifail, ir, is, it, itype, it1, it2, j, k, m, meshl
      integer n, nid, nit, nm, nrg, ptclsl
      logical azimth, resyn
      real rs, rf, strig, th
      include 'inc/pi.f'
c
      if( .not. p2d )call crash( 'P2GRNM', 'Wrong type of code' )
c open files
      if( ngt .lt. 0 )call new_unit( ngt )
      nit = ngt
      call opnfil( nit, 'grt', 'unformatted', 'old', 'seq', i )
      if( i .eq. 0 )then
c verify existing transformed file
        call grnhed( .true., .false., ifail )
        if( ifail .ne. 0 )call crash( 'P2GRNM',
     +                            'Wrong transformed Green fun file' )
        if( direct )then
          if( ngd .lt. 0 )call new_unit( ngd )
          nid = ngd
          call opnfil( nid, 'grd', 'unformatted', 'old', 'seq', i )
          if( i .eq. 1 )call crash( 'P2GRNM',
     +                              'Direct Green fn file not found' )
c verify existing direct file
          call grnhed( .true., .true., ifail )
          if( ifail .ne. 0 )call crash( 'P2GRNM',
     +                                 'Wrong direct Green fun file' )
        end if
        if( master )print *, ' Old files found and checked in P2GRNM'
      else if( master )then
c should not create this file on the slaves as well as on the master
        nrg = nr( jgrid )
c preserve value of lptcls
        ptclsl = lptcls
c allocate space
        meshl = 2 * nrg * ( na / 2 + 1 )
        allocate ( grdmss( meshl, 1 ) )
        lptcls = nrg * nrg * ng
        allocate ( ptcls( lptcls ) )
        lwork = 2 * nrg * ( na / 2 + 1 )
        allocate ( w( lwork ) )
        ltrig = 3 * na + 15
        allocate ( trig( ltrig ) )
c create new files
        print *, 'P2GRNM: Creating new Green fn files'
        call opnfil( nit, 'grt', 'unformatted', 'new', 'seq', i )
        call grnhed( .false., .false., ifail )
        if( direct )then
          if( ngd .lt. 0 )call new_unit( ngd )
          nid = ngd
          call opnfil( nid, 'grd', 'unformatted', 'new', 'seq', i )
          call grnhed( .false., .true., ifail )
        end if
c form radial forces, azimuthal forces and then potentials
        resyn = ng .lt. mmax
        do itype = 1, 3
          azimth = itype .eq. 2
c initialize FFTPAK routines
          if( azimth )then
            it1 = 2
            it2 = mmax - 1
            call rsinti( na / 2 - 1, trig )
            strig = 1. / real( na )
          else
            it1 = 1
            it2 = mmax
            call rcosti( na / 2 + 1, trig )
            strig = 1. / real( na )
          end if
c scan over radii of sources
          n = 1
          do is = 1, nrg
            m = 0
            rs = grofu( real( is - 1 ) )
c generate table of direct values
            do ir = 1, nrg
              j = ir
              rf = grofu( real( ir - 1 ) )
              do it = it1, it2
                th = 2. * pi * real( it - 1 ) / real( na )
                m = m + 1
                if( itype .eq. 1 )grdmss( m, 1 ) =
     +                                          radfor( rs, rf, th, 0. )
                if( itype .eq. 2 )grdmss( m, 1 ) =
     +                                          azifor( rs, rf, th, 0. )
                if( itype .eq. 3 )grdmss( m, 1 ) =
     +                                          grnfun( rs, rf, th, 0. )
                w( j ) = grdmss( m, 1 )
                j = j + nrg
              end do
            end do
            if( direct .and. ( .not. resyn ) )then
              write( nid )( grdmss( i, 1 ), i = 1, m )
            end if
            if( azimth )then
c find Fourier sine transform of azimuthal forces
              i = 1
              do ir = 1, nrg
                call rsint( na / 2 - 1, grdmss( i, 1 ), trig )
                i = i + na / 2 - 1
              end do
            else
c find Fourier cosine transform of radial forces or potentials
              i = 1
              do ir = 1, nrg
                call rcost( na / 2 + 1, grdmss( i, 1 ), trig )
                i = i + na / 2 + 1
              end do
            end if
c filter out unwanted harmonics
            m = 0
            do ir = 1, nrg
              j = ir
              do it = it1, it2
                m = m + 1
                if( it .gt. ng )grdmss( m, 1 ) = 0
                if( mod( it - 1, nsect ) .ne. 0 )grdmss( m, 1 ) = 0
                j = j + nrg
              end do
            end do
c store transformed matrix
            nm = ng
            k = 1 - mmax
            if( azimth )then
              nm = min( mmax - 2, ng - 1 )
              k = k + 2
            end if
            do i = 1, nrg
              k = k + mmax
              if( azimth )k = k - 2
              call blkcpy( grdmss( k, 1 ), ptcls( n ), nm )
              n = n + nm
            end do
c re-synthesise direct version
            if( direct .and. resyn )then
              if( azimth )then
                m = 1
                do ir = 1, nrg
                  call rsint( na / 2 - 1, grdmss( m, 1 ), trig )
                  m = m + na / 2 - 1
                end do
              else
c find Fourier cosine transform of radial forces or potentials
                m = 1
                do ir = 1, nrg
                  call rcost( na / 2 + 1, grdmss( m, 1 ), trig )
                  m = m + na / 2 + 1
                end do
              end if
c normalize and save new table of direct values
              m = m - 1
              do i = 1, m
                grdmss( i, 1 ) = strig * grdmss( i, 1 )
              end do
              write( nid )( grdmss( i, 1 ), i = 1, m )
            end if
          end do
c reformat transformed version
          call p2gwrt( azimth, lwork, w )
        end do
        deallocate ( ptcls )
        deallocate ( w )
        deallocate ( trig )
        deallocate ( grdmss )
c restore value of lptcls
        lptcls = ptclsl
c close and re-open file(s) to ensure that writes are completed
        close( nit )
        call opnfil( nit, 'grt', 'unformatted', 'old', 'seq', i )
        if( direct )then
          close( nid )
          call opnfil( nid, 'grd', 'unformatted', 'old', 'seq', i )
        end if
        write( no, '( '' Green function creation complete'' )' )
        print *, 'P2GRNM: Green fn files ready'
      end if
      return
      end
