      subroutine c2chkf
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to check the solution for the field on the 2-D Cartesian grid
c   by direct convolution - called from C2TEST
c an empty mass array except for a small number of unit masses should
c   have been first created by a call to C2CKST and the field for that
c   determined in the usual way by a call to C2FNDF
c this routine is the only routine to make use of the direct Green
c   function file created by C2GRNM
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
c
      real, allocatable :: w(:)
c
c local variables
      character*10 a( 3 )
      integer g, i, itype, j, k, kt, l, m, n, p, q, s
      real err, errmax, x, y
c
      data a / '  x forces', '  y forces', 'Potentials' /
c
c check Green function file header
      call grnhed( .true., .true., i )
      if( i .ne. 0 )call crash( 'C2CHKF',
     +                          'Incorrect direct Green function' )
c allocate scratch areas
      allocate ( w( mesh( jgrid ) ) )
      l = ( ngx + 1 ) * ( ngy + 1 )
      allocate ( grd( l ) )
c
      if( master )write( no, 200 )
c work over field types
      kt = 3
      if( .not. phys )kt = 2
      do itype = 1, kt
c clear results space
        do i = 1, mesh( jgrid )
          w( i ) = 0.
        end do
        if( itype .eq. 1 )then
          l = ngx * ( ngy + 1 )
        else if( itype .eq. 2 )then
          l = ( ngx + 1 ) * ngy
        else if( itype .eq. 3 )then
          l = ( ngx + 1 ) * ( ngy + 1 )
        end if
        read( ngd )( grd( i ), i = 1, l )
c convolve with all sources
        do s = 1, nmass
          l = r( s )
          m = th( s )
          g = 0
c work through Green function - x-accelerations
          if( itype .eq. 1 )then
            do j = 0, ngy
              k = ( m - 1 + j ) * ngx
              n = ( m - 1 - j ) * ngx
              do i = 1, ngx
                g = g + 1
                p = l + i
                q = l - i
                if( k .lt. mesh( jgrid ) )then
                  if( p .le. ngx )w( k + p ) = w( k + p ) + grd( g )
                  if( q .gt. 0 )w( k + q ) = w( k + q ) - grd( g )
                end if
                if( ( j .gt. 0 ) .and. ( n .ge. 0 ) )then
                  if( p .le. ngx )w( n + p ) = w( n + p ) + grd( g )
                  if( q .gt. 0 )w( n + q ) = w( n + q ) - grd( g )
                end if
              end do
            end do
c work through Green function - y-accelerations
          else if( itype .eq. 2 )then
            do j = 1, ngy
              k = ( m - 1 + j ) * ngx
              n = ( m - 1 - j ) * ngx
              do i = 0, ngx
                g = g + 1
                p = l + i
                q = l - i
                if( k .lt. mesh( jgrid ) )then
                  if( p .le. ngx )w( k + p ) = w( k + p ) + grd( g )
                  if( ( i .gt. 0 ) .and. ( q .gt. 0 ) )w( k + q ) =
     +                                         w( k + q ) + grd( g )
                end if
                if( n .ge. 0 )then
                  if( p .le. ngx )w( n + p ) = w( n + p ) - grd( g )
                  if( ( i .gt. 0 ) .and. ( q .gt. 0 ) )w( n + q ) =
     +                                         w( n + q ) - grd( g )
                end if
              end do
            end do
c work through Green function - potentials
          else if( itype .eq. 3 )then
            do j = 0, ngy
              k = ( m - 1 + j ) * ngx
              n = ( m - 1 - j ) * ngx
              do i = 0, ngx
                g = g + 1
                p = l + i
                q = l - i
                if( k .lt. mesh( jgrid ) )then
                  if( p .le. ngx )w( k + p ) = w( k + p ) + grd( g )
                  if( ( i .gt. 0 ) .and. ( q .gt. 0 ) )w( k + q ) =
     +                                         w( k + q ) + grd( g )
                end if
                if( ( j .gt. 0 ) .and. ( n .ge. 0 ) )then
                  if( p .le. ngx )w( n + p ) = w( n + p ) + grd( g )
                  if( ( i .gt. 0 ) .and. ( q .gt. 0 ) )w( n + q ) =
     +                                         w( n + q ) + grd( g )
                end if
              end do
            end do
          end if
        end do
c comparison
        errmax = 0.
        k = -99
        l = -99
        m = 0
        do j = 1, ngy
          do i = 1, ngx
            m = m + 1
            if( abs( w( m ) ) .gt. .1e-7 )then
              err = abs( grdfld( m, itype ) - w( m ) )
              if( err .gt. errmax )then
                errmax = err
                x = grdfld( m, itype )
                y = w( m )
                k = i
                l = j
              end if
            end if
          end do
        end do
        if( master )write( no, 201 )a( itype ), errmax, k, l, x, y
      end do
      deallocate ( grd )
      deallocate ( w )
  200 format( /' Field type  Max abs error   x   y  values' )
  201 format( a10, e16.5, 2i4, 2f14.9 )
      return
      end
