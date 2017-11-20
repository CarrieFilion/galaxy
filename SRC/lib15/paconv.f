      subroutine paconv( itype, lw, mt, ft )
c  Copyright (C) 2015, Jerry Sellwood
c
c performs the radial convolution part of the solution for the fields
c   on the axisymmetric polar grid
c The Green function is presumed stored in Fourier transformed form
c   on a file (created by PAGRNM) and is read in to be used here.  If
c   the local array green is dimensioned to be large enough, it is
c   read in once at the first call and stored, otherwise it is read in
c   record-by-record as needed, at each call.
c The routine is called for each type of convolution, for radial forces,
c   vertical forces and potentials.  The transformed mass array is
c   assuemd to start at msh( ip( 4 ) ) and the results are copied at
c   the end to msh( ip( itype) ) - NB the potential coefficients will
c   overwrite the first part of the mass coefficients!
      use aarrays
      implicit none
c
c calling arguments
      integer itype, lw
      real mt( lw ), ft( lw )
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
      real, allocatable :: gr(:)
c
      real, allocatable, save :: green( :, :, : )
c
c local variables
      integer i, ifail, ir, is, iz, iuse, j, jty, jz, k, kty, kz, m
      integer nrg
      logical file, nskip
      save iuse, file
c
      data iuse / 0 /
c
      nrg = nr( jgrid )
c check Green function file header
      if( iuse .eq. 0 )then
        rewind ngt
        call grnhed( .true., .false., ifail )
        if( ifail .ne. 0 )call crash( 'PACONV',
     +                                     'Wrong Green function file' )
c file should be false when running in parallel to avoid collisions
        file = .not. parallel
        if( file )then
c check workspace to be used instead
          if( master )then
            print *, 'Green fn used from file'
            write( no, * )'Green fn used from file'
          end if
        else
c read in Green function
          if( master )then
            print *, 'Green function stored in memory'
            write( no, * )'Green function stored in memory'
          end if
          allocate( green( nrg * nrg, ngz + 1, 3 ) )
          if( fixrad )then
            jty = 2
            kty = 2
          else
            jty = 1
            kty = 3
          end if
          do j = jty, kty
            jz = 1
            kz = ngz + 1
            if( j .eq. 2 )then
              jz = 2
              kz = ngz
            end if
            k = nrg * nrg
            do iz = jz, kz
              read( ngt )( green( i, iz, j ), i = 1, k )
            end do
          end do
          close( ngt )
        end if
        iuse = 1
      end if
c allocate space
      j = nrg * nrg
      allocate ( gr( j ) )
c
c skip header if necessary
      if( file .and. ( ( itype .eq. 1 ) .or.
     +               ( ( itype .eq. 2 ) .and. fixrad ) ) )then
        rewind ngt
        read( ngt )
      end if
c
      jz = 1
      kz = ngz + 1
      if( itype .eq. 2 )then
        jz = 2
        kz = ngz
      end if
c clear results array
      do i = 1, 2 * mesh( jgrid )
        ft( i ) = 0
      end do
c convolve with each Fourier component in turn
      do iz = jz, kz
        nskip = ( iz .gt. 1 ) .and. ( iz .le. ngz )
        if( file )then
          k = nrg * nrg
          read( ngt )( gr( i ), i = 1, k )
        end if
        j = 0
c sum over radii for sources
        do is = 1, nrg
c k & m point to the cosine terms, the sine terms are k + 1 and m + 1
          k = max( 1, 2 * ( iz - 1 ) )
          m = 2 * ngz * ( is - 1 ) + k
          if( itype .ne. 2 )then
c work over radial field points - radial forces or potentials are symmetric
            if( file )then
              do ir = 1, nrg
                j = j + 1
                ft( k ) = ft( k ) + mt( m ) * gr( j )
                if( nskip )ft( k + 1 ) =
     +                               ft( k + 1 ) + mt( m + 1 ) * gr( j )
                k = k + 2 * ngz
              end do
            else
              do ir = 1, nrg
                j = j + 1
                ft( k ) = ft( k ) + mt( m ) * green( j, iz, itype )
                if( nskip )ft( k + 1 ) = ft( k + 1 ) +
     +                         mt( m + 1 ) * green( j, iz, itype )
                k = k + 2 * ngz
              end do
            end if
          else
c work over radial field points - vertical forces are anti-symmetric
            if( file )then
              do ir = 1, nrg
                j = j + 1
                ft( k ) = ft( k ) + mt( m + 1 ) * gr( j )
                ft( k + 1 ) = ft( k + 1 ) - mt( m ) * gr( j )
                k = k + 2 * ngz
              end do
            else
              do ir = 1, nrg
                j = j + 1
                ft( k ) = ft( k ) + mt( m + 1 ) * green( j, iz, itype )
                ft( k + 1 ) = ft( k + 1 ) -
     +                                  mt( m ) * green( j, iz, itype )
                k = k + 2 * ngz
              end do
            end if
          end if
c end loop over sources
        end do
c end loop over Fourier components
      end do
c return local allocated space
      deallocate ( gr )
      return
      end
