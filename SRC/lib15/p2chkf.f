      subroutine p2chkf
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to check the solution for the field on the 2-D polar code
c   by direct convolution - called from P2TEST
c an empty mass array except for a small number of unit masses should
c   have been first created by a call to P2CKST and the field for that
c   determined in the usual way by a call to P2FNDF
c this routine is the only routine to make use of the direct Green
c   function file created by P2GRNM
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
c local allocatable arrays
      real, allocatable :: grd(:)
      real, allocatable :: w(:)
c
c local variables
c      character*12 a( 3 )
      integer i, ir, it, itype, j, jrf, k, krf, kt, l, m, mr, mt, n, nt
      logical azimth
      real err, errmax, x, y
c
c      data a / 'Radial forcs', 'Azimuthal fs', 'Potentials  ' /
c
c check Green function file header
      call grnhed( .true., .true., i )
      if( i .ne. 0 )call crash( 'P2CHKF',
     +                          'Incorrect direct Green function' )
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
      write( no, 200 )
c allocate space
      j = nr( jgrid ) * mmax
      allocate ( grd( j ) )
c allocate only active field region
      j = ( jrf - 1 ) * na + 1
      k = krf * na
      allocate ( w( j:k ) )
c work over field types
      kt = 3
      if( .not. phys )kt = 2
      itype = 0
      if( fixrad )itype = 1
      do while ( itype .lt. kt )
        itype = itype + 1
c azimuthal forces are all zero if ng = 1
        if( ( itype .eq. 2 ) .and. ( ng .eq. 1 ) )itype = itype + 1
        azimth = itype .eq. 2
        nt = mmax
        if( azimth )nt = mmax - 2
c
c clear results space
        j = ( jrf - 1 ) * na + 1
        k = krf * na
        do i = j, k
          w( i ) = 0.
        end do
c work through grid rings for sources
        mt = mmax - 2
        do mr = 1, nr( jgrid )
c read in direct Green function
          k = nr( jgrid ) * nt
          read( ngd )( grd( i ), i = 1, k )
c convolve with any masses on this ring
          do n = 1, nmass
            if( r( n ) .eq. mr )then
              m = th( n )
              l = ( jrf - 1 ) * na
              if( azimth )then
                j = ( jrf - 1 ) * mt
              else
                j = ( jrf - 1 ) * mmax + 1
              end if
c sum over radii for this ring
              do ir = jrf, krf
c work over angles
                k = l + m
                if( .not. azimth )w( k ) = w( k ) + grd( j )
                do it = 1, mt
                  j = j + 1
                  k = m + it
                  if( k .gt. na )k = k - na
                  k = k + l
                  w( k ) = w( k ) + grd( j )
                  k = m - it
                  if( k .le. 0 )k = k + na
                  k = k + l
                  if( azimth )then
                    w( k ) = w( k ) - grd( j )
                  else
                    w( k ) = w( k ) + grd( j )
                  end if
                end do
                if( .not. azimth )then
                  k = m + mmax - 1
                  if( k .gt. na )k = k - na
                  k = k + l
                  w( k ) = w( k ) + grd( j + 1 )
                  j = j + 2
                end if
                l = l + na
              end do
            end if
          end do
        end do
c comparison
        errmax = 0.
        k = -99
        l = -99
        n = ( jrf - 1 ) * na
        do ir = jrf, krf
          do it = 1, na
            n = n + 1
            if( abs( w( i ) ) .gt. .1e-7 )then
              err = abs( grdfld( n, itype ) - w( n ) )
              if( err .gt. errmax )then
                errmax = err
                x = grdfld( n, itype )
                y = w( n )
                k = ir
                l = it
              end if
            end if
          end do
        end do
        write( no, 201 )itype, errmax, k, l, x, y
c end loop over field types
      end do
      return
  200 format( /' Field type  Max abs error   r   t  values' )
  201 format( i10, e16.5, 2i4, 2f14.9 )
      end
