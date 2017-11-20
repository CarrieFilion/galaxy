      subroutine c2conv( itype, l, msstr )
c  Copyright (C) 2015, Jerry Sellwood
c
c multiplies the transformed mass array by the appropriate Green function
c   coefficients in the solution for the field on the 2-D Cartesian grid
c The Green function is presumed stored in Fourier transformed form
c   in a file (created by C2GRNM) and is read in to be used here.  If
c   the local array green is dimensioned to be large enough, it is
c   read in once at the first call and stored, otherwise it is read in
c   record-by-record as needed, at each call.
c The routine is called for each type of convolution, for x-accelerations,
c   y-accelerations and potentials
c The results over-write the input array
      use aarrays
      implicit none
c
c calling arguments
      integer itype, l
      real msstr( l )
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
      real, allocatable :: gr(:)
c
c local variables
      integer i, ifail, ix, iy, iuse, j, jx, jy, k, kx, ky
      logical nyskip
      real x
      save iuse
c
      data iuse / 0 /
c
c check Green function file header
      if( iuse .eq. 0 )then
        call grnhed( .true., .false., ifail )
        if( ifail .ne. 0 )call crash( 'C2CONV',
     +                                     'Wrong Green function file' )
        if( master )then
          print *, 'Green fn used from file'
          write( no, * )'Green fn used from file'
        end if
        iuse = 1
      end if
      j = ngx + 1
      allocate ( gr( j ) )
c skip header if necessary
      if( itype .eq. 1 )then
        rewind ngt
        read( ngt )
      end if
c set constants
      if( itype .eq. 1 )then
        jx = 2
        kx = ngx
      else
        jx = 1
        kx = ngx + 1
      end if
      if( itype .eq. 2 )then
        jy = 2
        ky = ngy
      else
        jy = 1
        ky = ngy + 1
      end if
c zero implied for end terms of sine expansion of the Green function
      if( itype .eq. 1 )then
        k = 0
        do iy = 1, 2 * ngy
          k = k + 1
          msstr( k ) = 0
          k = k + 2 * ngx - 1
          msstr( k ) = 0
        end do
      else if( itype .eq. 2 )then
        k = 0
        do ix = 1, 2 * ngx
          k = k + 1
          msstr( k ) = 0
        end do
        k = 2 * ngx * ( 2 * ngy - 1 )
        do ix = 1, 2 * ngx
          k = k + 1
          msstr( k ) = 0
        end do
      end if
c work over y Fourier components
      do iy = jy, ky
c get Green function
        read( ngt )( gr( ix ), ix = jx, kx )
c i & k point to the cc & cs terms respectively, i + 1 & k + 1 to sc & ss
        i = max( 1, 2 * ( iy - 1 ) )
        i = ( i - 1 ) * 2 * ngx + 1
        k = i + 2 * ngx
c work over x Fourier components of x-accelerations
        if( itype .eq. 1 )then
          nyskip = ( iy .gt. 1 ) .and. ( iy .le. ngy )
          i = i + 1
          k = k + 1
          do ix = 2, ngx
            x = msstr( i + 1 ) * gr( ix )
            msstr( i + 1 ) = -msstr( i ) * gr( ix )
            msstr( i     ) =  x
            if( nyskip )then
              x = msstr( k + 1 ) * gr( ix )
              msstr( k + 1 ) = -msstr( k ) * gr( ix )
              msstr( k     ) =  x
            end if
            i = i + 2
            k = k + 2
          end do
c work over x Fourier components of y-accelerations
        else if( itype .eq. 2 )then
          do ix = 1, ngx + 1
            x = msstr( k ) * gr( ix )
            msstr( k ) = -msstr( i ) * gr( ix )
            msstr( i ) =  x
            if( ( ix .gt. 1 ) .and. ( ix .le. ngx ) )then
              x = msstr( k + 1 ) * gr( ix )
              msstr( k + 1 ) = -msstr( i + 1 ) * gr( ix )
              msstr( i + 1 ) =  x
              i = i + 2
              k = k + 2
            else
              i = i + 1
              k = k + 1
            end if
          end do
c work over x Fourier components of potentials
        else if( itype .eq. 3 )then
          nyskip = ( iy .gt. 1 ) .and. ( iy .le. ngy )
          do ix = 1, ngx + 1
            msstr( i ) = msstr( i ) * gr( ix )
            if( nyskip )msstr( k ) = msstr( k ) * gr( ix )
            if( ( ix .gt. 1 ) .and. ( ix .le. ngx ) )then
              msstr( i + 1 ) = msstr( i + 1 ) * gr( ix )
              if( nyskip )msstr( k + 1 ) = msstr( k + 1 ) * gr( ix )
              i = i + 2
              k = k + 2
            else
              i = i + 1
              k = k + 1
            end if
          end do
        end if
      end do
      deallocate ( gr )
      return
      end
