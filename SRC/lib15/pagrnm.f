      subroutine pagrnm( direct )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c sets up Green function file for 3-D axisymmetric grid or checks that
c   a pre-existing file is appropriate for the current grid
c   uses QUADPACK routine DQAWC and Burkardt's routine QUAD
c
c this version uses single precision FFTPAK routines
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
      common / feildp / rf, zf
      real*8 rf, zf
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c externals
      external pafrfc, pafrfn, pafzfn, papotf, papotc
      real grofu
      real*8 quad_chy, quad_Pat
c
c local arrays
      real, allocatable :: trig( : )
      real, allocatable :: wrk( : )
      real, allocatable :: x( : )
c
c local variables
      integer i, ier, is, ir, iz, itype, j, jty, k, kty, m, n
      integer nid, nit, nrg, nz
      real*8 ar, az, epsa, epsr, phi, rs, r1, r2
c
      nrg = nr( jgrid )
      if( ngt .lt. 0 )call new_unit( ngt )
      nit = ngt
      if( direct )then
c open file
        if( ngd .lt. 0 )call new_unit( ngd )
        nid = ngd
        call opnfil( nid, 'grd', 'unformatted', 'old', 'seq', i )
        if( i .eq. 0 )then
c verify existing file
          call grnhed( .true., .true., ier )
          if( ier .ne. 0 )call crash( 'PAGRNM',
     +                              'Wrong direct Green function file' )
          print *, ' Old direct file found and checked in PAGRNM'
        end if
      end if
c open file
      call opnfil( nit, 'grt', 'unformatted', 'old', 'seq', i )
      if( i .eq. 0 )then
c verify existing file
        call grnhed( .true., .false., ier )
        if( ier .ne. 0 )call crash( 'PAGRNM',
     +                           'Wrong transformed Green functn file' )
        print *, ' Old transformed file found and checked in PAGRNM'
      else if( master )then
c should not create this file on the slaves as well as on the master
        print *, 'PAGRNM: Creating new Green fn file'
c allocate space
        allocate ( x( ( ngz + 1 ) * ( ngz + 1 ) ) )
        allocate ( trig( 3 * ( ngz + 1 ) + 15 ) )
        n = nrg * nrg * ( ngz + 1 )
        allocate ( wrk( n ) )
c open Green function file(s) and write header(s)
        call opnfil( nit, 'grt', 'unformatted', 'new', 'seq', i )
        call grnhed( .false., .false., ier )
        if( direct )then
          call opnfil( nid, 'grd', 'unformatted', 'new', 'seq', i )
          call grnhed( .false., .true., ier )
        end if
        if( fixrad )then
          jty = 2
          kty = 2
        else
          jty = 1
          kty = 3
        end if
        do itype = jty, kty
          nz = ngz + 1
          if( itype .eq. 2 )nz = ngz
          m = 0
c generate tables of values
          do iz = 1, nz
            zf = dzg * real( iz - 1 )
            do is = 1, nrg
              rs = grofu( real( is - 1 ) )
              do ir = 1, nrg
                rf = grofu( real( ir - 1 ) )
c compute attraction of a uniform surface density annulus which extends
c   half a mesh space on either side
                r1 = max( real( is ) - 1.5, 0. )
                r1 = grofu( sngl( r1 ) )
                r2 = grofu( real( is ) - 0.5 )
c radial forces
                if( itype .eq. 1 )then
c need a Cauchy principal value integrator for radial part of self-force
                  if( ( is .gt. 1 ) .and. ( ir .eq. is ) .and.
     +                ( iz .eq. 1 ) )then
                    epsa = 1.d-5
                    epsr = epsa
                    ar = quad_chy( pafrfc, r1, r2, rs, epsa, epsr, ier )
                    if( ier .gt. 0 )then
                      print *, 'ier =', ier, ' from QUAD_CHY'
                      call crash( 'PAGRNM', 'QUADPACK error' )
                    end if
                  else
c integrand is regular
                    ar = quad_Pat( pafrfn, r1, r2, epsa, epsr, ier )
                  end if
                  ar = 2. * ar / ( r2**2 - r1**2 )
                end if
c vertical forces
                if( itype .eq. 2 )then
                  if( iz .eq. 1 )then
                    az = 0
                  else
                    if( ir .ge. is )then
c integrand is regular everywhere
                      epsa = 1.d-8
                      epsr = epsa
                      az = quad_Pat( pafzfn, r1, r2, epsa, epsr, ier )
                      az = 2. * az / ( r2**2 - r1**2 )
                    else
c impose equal and opposite vertical forces
                      n = nrg * ( nrg * ( iz - 1 ) + ir - 1 ) + is
                      az = wrk( n )
                    end if
                  end if
                end if
c potential
                if( itype .eq. 3 )then
c need a Cauchy principal value integrator for potential
                  if( ( is .gt. 1 ) .and. ( ir .eq. is ) .and.
     +                ( iz .eq. 1 ) )then
                    epsa = 1.d-5
                    epsr = epsa
                    phi =
     +                   quad_chy( papotc, r1, r2, rs, epsa, epsr, ier )
                    if( ier .gt. 0 )then
                      print *, 'ier =', ier, ' from QUAD_CHY'
                      call crash( 'PAGRNM', 'QUADPACK error' )
                    end if
                  else
c integrand is regular
                    epsa = 1.d-8
                    epsr = epsa
                    phi = quad_Pat( papotc, r1, r2, epsa, epsr, ier )
                  end if
                  phi = 2. * phi / ( r2**2 - r1**2 )
                end if
c store these values
                m = m + 1
                if( itype .eq. 1 )wrk( m ) = ar
                if( itype .eq. 2 )wrk( m ) = az
                if( itype .eq. 3 )wrk( m ) = phi
              end do
            end do
          end do
          if( itype .eq. 1 )print *, 'Created table of radial forces'
          if( itype .eq. 2 )print *, 'Created table of vertical forces'
          if( itype .eq. 3 )print *, 'Created table of potentials'
c save direct values if requested
          m = nrg * nrg
          if( direct )then
            k = 0
            do iz = 1, ngz
              j = k + 1
              k = k + m
              write( nid )( wrk( i ), i = j, k )
            end do
          end if
c initialize FFTs
          if( itype .eq. 2 )then
            call rsinti( ngz - 1, trig )
            nz = ngz - 1
          else
            call rcosti( ngz + 1, trig )
            nz = ngz + 1
          end if
c copy, transform and restore columns of values - it would require too
c  much memory to do this in one operation per itype
          k = 0
          do ir = 1, m
            j = k + ir
            if( itype .eq. 2 )j = j + m
            do iz = 1, nz
              x( iz ) = wrk( j )
              j = j + m
            end do
            if( itype .eq. 2 )then
              call rsint( ngz - 1, x, trig )
            else
              call rcost( ngz + 1, x, trig )
            end if
c save transformed values
            j = k + ir
            if( itype .eq. 2 )j = j + m
            do iz = 1, nz
              wrk( j ) = x( iz )
              j = j + m
            end do
          end do
c write out transformed function
          k = 0
          if( itype .eq. 2 )k = k + m
          do iz = 1, nz
            j = k + 1
            k = k + m
            write( nit )( wrk( i ), i = j, k )
          end do
          if( master )print *, 'Fourier coefficients written to file'
        end do
c close and re-open file(s) to ensure that writes are completed
        close( nit )
        call opnfil( nit, 'grt', 'unformatted', 'old', 'seq', i )
        if( direct )then
          close( nid )
          call opnfil( nid, 'grd', 'unformatted', 'old', 'seq', i )
        end if
        print *, 'Green function creation complete'
        write( no, '( '' Green function creation complete'' )' )
        print *, 'PAGRNM: Green fn file ready'
      end if
      return
      end
