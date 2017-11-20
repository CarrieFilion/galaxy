      subroutine pachkf
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to perform a cross check of the FFT solution for the
c   field of an arbitrary mass array on the axisymmetric polar grid
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
      real, allocatable :: green(:,:,:)
c
      real, allocatable :: g(:)
c
      real, allocatable :: wl(:,:)
c
c local variables
      character*12 a( 3 )
      logical file
      integer i, ir, itype, iz, j, jr, jz, jty, k, kty, kz, l, n
      integer nrg
      real err, errmax, sign, x, y
c
      data a / 'Radial forcs', 'Vertcl forcs', 'Potentials  ' /
c
c allocate space
      nrg = nr( jgrid )
      allocate ( wl( nrg * ngz, 3 ) )
c
      if( fixrad )then
        jty = 2
        kty = 2
      else
        jty = 1
        kty = 3
      end if
c save current results
      if( fixrad )then
        jty = 2
        kty = 2
      else
        jty = 1
        kty = 2
        if( phys )kty = 3
      end if
c check Green function file header
      rewind ngd
      call grnhed( .true., .true., i )
      if( i .ne. 0 )call crash( 'PACHKF',
     +                            'Wrong Green function file' )
c determine whether green array is large enough
      file = .false.
      if( file )then
c check workspace to be used instead
        print *, 'Green fn used from file'
        allocate ( g( nrg * nrg ) )
      else
c read in Green function
        print *, 'Green function stored in memory'
        allocate ( green( nrg * nrg, ngz, 3 ) )
        do j = jty, kty
          do iz = 1, ngz
            read( ngd )( green( i, iz, j ), i = 1, nrg * nrg )
          end do
        end do
        close( ngd )
      end if
c skip header if necessary
      if( file )then
        rewind ngd
        read( ngd )
      end if
c work over types
      do itype = jty, kty
c clear results area
        do j = 1, mesh( jgrid )
          wl( j, itype ) = 0
        end do
c work over planes of Green's function
        do jz = 1, ngz
          if( file )then
            k = nrg * nrg
            read( ngd )( g( i ), i = 1, k )
          end if
c work over sources
          do i = 1, nmass
c planes above and containing the source
            k = ( r( i ) - 1 ) * nrg
            kz = z( i ) + jz - 2
            if( kz .lt. ngz )then
              l = kz * nrg
              if( file )then
                do jr = 1, nrg
                  k = k + 1
                  l = l + 1
                  wl( l, itype ) = wl( l, itype ) + g( k )
                end do
              else
                do jr = 1, nrg
                  k = k + 1
                  l = l + 1
                  wl( l, itype ) =
     +                            wl( l, itype ) + green( k, jz, itype )
                end do
              end if
            end if
c planes below the source
            sign = 1
            if( itype .eq. 2 )sign = -1
            kz = z( i ) - jz
            if( ( jz. gt. 1 ) .and. ( kz .ge. 0 ) )then
              k = ( r( i ) - 1 ) * nrg
              l = kz * nrg
              if( file )then
                do jr = 1, nrg
                  k = k + 1
                  l = l + 1
                  wl( l, itype ) = wl( l, itype ) + sign * g( k )
                end do
              else
                do jr = 1, nrg
                  k = k + 1
                  l = l + 1
                  wl( l, itype ) =
     +                     wl( l, itype ) + sign * green( k, jz, itype )
                end do
              end if
            end if
c end sum over sources
          end do
c end loop over planes
        end do
c end loop over types
      end do
c
c comparison
      if( master )write( no, * )
      if( master )write( no, * )'Results from PACHKF:'
c work over types
      do itype = jty, kty
        n = 0
        errmax = 0.
        k = -99
        l = -99
        do iz = 1, ngz
          do ir = 1, nr( jgrid )
            n = n + 1
            err = abs( grdfld( n, itype ) - wl( n, itype ) )
            if( err .gt. errmax )then
              errmax = err
              x = wl( n, itype )
              y = grdfld( n, itype )
              k = ir
              l = iz
            end if
          end do
        end do
        if( master )write( no, '( / 5x, a12 )' )a( itype )
        if( master )write( no, 200 )errmax, k, l, x, y
      end do
c return local workspace
      if( file )then
        deallocate ( g )
      else
        deallocate ( green )
      end if
      deallocate ( wl )
      return
  200 format( ' Max absolute error', e12.4, ' at', 2i4,
     +        ' for values', 2f12.8 )
      end
