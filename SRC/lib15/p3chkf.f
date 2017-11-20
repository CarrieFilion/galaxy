      subroutine p3chkf
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to check the solution for the field on the 3-D polar code by
c   direct convolution for a few point masses.  It is called from P3TEST
c An empty mass array except for a small number of unit masses should
c   have been first created by a call to P3CKST and the field for that
c   determined in the usual way by a call to P3FNDF
c This routine is the only routine to make use of the direct Green
c   function file created by P3GRNM
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
      include 'inc/source.f'
c
c local allocatable array
      real, allocatable :: grd(:)
      real, allocatable :: w(:)
c
c local variables
c      character*12 a( 4 )
      integer i, ir, is, it, itype, iz, iz1, iz2, j, jrf, k
      integer krf, kt, l, m, n, nt, nz
      logical azimth, linc, lower, upper, vert
      real err, errmax, sgna, sgnv, x, y
c
c      data a / 'Radial forcs', 'Azmthal frcs',
c     +         'Vertical frc', 'Potentials  ' /
c
c check Green function file header
      call grnhed( .true., .true., i )
      if( i .ne. 0 )call crash( 'P3CHKF',
     +                          'Incorrect direct Green function' )
c
c same condition as in findf
      if( potl .or. ( mzone .eq. nzones ) .or. wholeg )then
        jrf = 1
        krf = nr( jgrid )
      else
        jrf = nr( jgrid )
        krf = 1
        do i = 1, nzones
          jrf = min( jrf, jrad( i ) )
          krf = max( krf, krad( i ) )
        end do
      end if
      if( master )write( no, 200 )
c allocate space
      j = nr( jgrid ) * ( mmax + 1 ) * ( ngz + 1 )
      allocate ( grd( j ) )
c allocate only active field region
      j = ( jrf - 1 ) * na * ngz + 1
      k = krf * na * ngz
      allocate ( w( j:k ) )
c work over field types
      kt = 4
      if( .not. phys )kt = 3
      itype = 0
      do while ( itype .lt. kt )
        itype = itype + 1
c azimuthal forces are all zero if ng = 1
        if( ( itype .eq. 2 ) .and. ( ng .eq. 1 ) )itype = itype + 1
        azimth = itype .eq. 2
        vert = itype .eq. 3
        nt = mmax
        if( azimth )nt = mmax - 2
        iz1 = 1
        iz2 = ngz + 1
        if( vert )then
          iz1 = 2
          iz2 = ngz
        end if
        nz = iz2 - iz1 + 1
        sgna = 1
        if( azimth )sgna = -1
        sgnv = 1
        if( vert )sgnv = -1
c
c clear results space
        j = ( jrf - 1 ) * na * ngz + 1
        k = krf * na * ngz
        do i = j, k
          w( i ) = 0
        end do
c work through mesh for sources
        do is = 1, nr( jgrid )
c read in direct Green function
          k = nr( jgrid ) * nt * nz
          read( ngd )( grd( i ), i = 1, k )
c convolve with any masses on this ring
          do n = 1, nmass
            if( r( n ) .eq. is )then
              j = 0
c work over field radii for this source
              do ir = 1, nr( jgrid )
                linc = ( ir .ge. jrf ) .and. ( ir .le. krf )
c work over grid planes
                do iz = iz1, iz2
                  l = z( n ) + iz - 2
                  k = th( n )
                  if( .not. azimth )j = j + 1
                  upper = l .lt. ngz
                  if( upper )then
                    l = na * ( l + ( ir - 1 ) * ngz )
                    if( linc .and. ( .not. azimth )
     +                               )w( k + l ) = w( k + l ) + grd( j )
                  end if
                  m = z( n ) - iz
                  lower = ( m .ge. 0 ) .and. ( iz .gt. 1 )
                  if( lower )then
                    m = na * ( m + ( ir - 1 ) * ngz )
                    if( linc .and. ( .not. azimth )
     +                        )w( k + m ) = w( k + m ) + sgnv * grd( j )
                  end if
c work over angles
                  do it = 2, mmax - 1
                    j = j + 1
                    k = th( n ) + it - 1
                    if( k .gt. na )k = k - na
                    if( linc )then
                      if( upper )w( k + l ) = w( k + l ) + grd( j )
                      if( lower )w( k + m ) = w( k + m ) +
     +                                              sgnv * grd( j )
                    end if
                    k = th( n ) + 1 - it
                    if( k .le. 0 )k = k + na
                    if( linc )then
                      if( upper )w( k + l ) = w( k + l ) +
     +                                                   sgna * grd( j )
                      if( lower )w( k + m ) = w( k + m ) +
     +                                            sgna * sgnv * grd( j )
                    end if
c end loop over angles
                  end do
                  if( .not. azimth )then
                    k = th( n ) + mmax - 1
                    if( k .gt. na )k = k - na
                    j = j + 1
                    if( linc )then
                      if( upper )w( k + l ) = w( k + l ) + grd( j )
                      if( lower )w( k + m ) = w( k + m ) +
     +                                              sgnv * grd( j )
                    end if
                  end if
c end loop over planes
                end do
c end loop over field radii
              end do
c end conditional loop over masses
            end if
          end do
c end loop over source radii
        end do
c comparison
        errmax = 0.
        k = -99
        l = -99
        m = -99
        i = ( jrf - 1 ) * na * ngz
        do ir = jrf, krf
          do iz = 1, ngz
            do it = 1, na
              i = i + 1
              if( abs( w( i ) ) .gt. .1e-7 )then
                err = abs( grdfld( i, itype ) - w( i ) )
                if( err .gt. errmax )then
                  errmax = err
                  x = grdfld( i, itype )
                  y = w( i )
                  l = ir
                  m = iz
                  k = it
                end if
              end if
            end do
          end do
        end do
        if( master )write( no, 201 )itype, errmax, l, m, k, x, y
c end loop over field types
      end do
      return
  200 format( /' Field type  Max abs error   r   z   t  values' )
  201 format( i10, e16.5, 3i4, 2f14.9 )
      end
