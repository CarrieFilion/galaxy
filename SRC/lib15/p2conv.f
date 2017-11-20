      subroutine p2conv( itype, jrs, krs, jrf, krf )
c  Copyright (C) 2015, Jerry Sellwood
c
c performs the radial convolution part of the solution for the fields
c   on the 2-D polar grid
c The Green function is presumed stored in Fourier transformed form
c   on a file (created by P2GRNM) and is read in to be used here.  If
c   the local array green is dimensioned to be large enough, it is
c   read in once at the first call and stored, otherwise it is read in
c   record-by-record as needed, at each call.
c The routine is called for each type of convolution, for radial forces,
c   azimuthal forces and potentials and the results overwrite the input
c   data array.  The parameters jrs & krs mark in the inner and outer
c   radial rings where masses have been defined (ie a partial grid) and
c   jrf, krf mark the range of rings for which the solution is required
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
c local allocatable arrays
      integer lgr
      real, allocatable :: gr(:)
c
      real, allocatable :: green( :, :, : )
      save green
c
      integer lwork
      real, allocatable :: w(:)
c
c local variables
      integer i, im, im1, im2, ips, ir, is, iuse, j, k, l, m, n, nrg
      logical azimth, file, nskip
      save file, iuse
c
      data lg / lnfilt * .false. /, iuse / 0 /
c
      if( ( itype .lt. 1 ) .or. ( itype .gt. 3 ) )then
        if( master )print *, 'itype = ', itype
        call crash( 'P2CONV', 'Impossible itype' )
      end if
      nrg = nr( jgrid )
c first call
      if( iuse .eq. 0 )then
c check Green function file header
        call grnhed( .true., .false., i )
        if( i .ne. 0 )call crash( 'P2CONV',
     +                            'Incorrect Green function file' )
c determine whether green array is large enough
        file = .true.
c        file = .false.
        if( .not. file )then
          if( master )then
            print *, 'Green function stored in memory'
            write( no, * )'Green function stored in memory'
          end if
c allocate permanent buffer for full Green fn
          allocate ( green( nrg * nrg, ng, 3 ) )
c read in Green function
          do j = 1, 3
            azimth = j .eq. 2
            do im = 1, ng
              nskip = ( im .gt. 1 ) .and. ( im .lt. mmax )
              if( ( .not. azimth ) .or. nskip )then
                read( ngt )( green( i, im, j ), i = 1, nrg * nrg )
              end if
            end do
          end do
          close( ngt )
        else
          if( master )then
            print *, 'Green fn used from file'
            write( no, * )'Green fn used from file'
          end if
        end if
        iuse = 1
      end if
c allocate local buffer for Green fn and work space
      lgr = nrg * nrg
      allocate ( gr( lgr ) )
      lwork = mesh( jgrid )
      allocate ( w( lwork ) )
c skip header if necessary
      if( ( itype .eq. 1 ) .and. file )then
        rewind ngt
        read( ngt )
      end if
c clear results area
      k = ( jrf - 1 ) * na + 1
      l = krf * na
      do i = k, l
        w( i ) = 0
      end do
      ips = na * ( jrs - 1 )
c work over over Fourier components
      azimth = itype .eq. 2
      im1 = 1
      im2 = ng
      if( azimth )then
        im1 = 2
        im2 = min( mmax - 1, ng )
      end if
      do im = im1, im2
c read in Green function
        if( file )then
          k = nrg * nrg
          read( ngt )( gr( i ), i = 1, k )
        end if
c apply Fourier filter
        if( .not. lg( im ) )then
          nskip = ( im .gt. 1 ) .and. ( im .lt. mmax )
c work over radial field points
          do ir = jrf, krf
c set pointers - k & l are for cos & sin terms of field
            k = na * ( ir - 1 ) + max( 1, 2 * ( im - 1 ) )
            l = k + 1
c m & n are for cos & sin terms of mass transform
            m = ips + max( 1, 2 * ( im - 1 ) )
            n = m + 1
c j is index for Green fn
            j = nrg * ( ir - 1 ) + jrs - 1
            if( .not. azimth )then
c sum over sources - radial forces or potentials
              if( file )then
                do is = jrs, krs
                  j = j + 1
                  w( k ) = w( k ) + gr( j ) * grdfld( m, itype )
                  if( nskip )w( l ) =
     +                     w( l ) + gr( j ) * grdfld( n, itype )
                  m = m + na
                  n = n + na
                end do
              else
                do is = jrs, krs
                  j = j + 1
                  w( k ) =
     +               w( k ) + green( j, im, itype ) * grdfld( m, itype )
                  if( nskip )w( l ) =
     +               w( l ) + green( j, im, itype ) * grdfld( n, itype )
                  m = m + na
                  n = n + na
                end do
              end if
            else
c sum over sources - azimuthal forces
              if( file )then
                do is = jrs, krs
                  j = j + 1
                  w( k ) = w( k ) + gr( j ) * grdfld( n, itype )
                  w( l ) = w( l ) - gr( j ) * grdfld( m, itype )
                  m = m + na
                  n = n + na
                end do
              else
                do is = jrs, krs
                  j = j + 1
                  w( k ) =
     +               w( k ) + green( j, im, itype ) * grdfld( n, itype )
                  w( l ) =
     +               w( l ) - green( j, im, itype ) * grdfld( m, itype )
                  m = m + na
                  n = n + na
                end do
              end if
            end if
c end loop over source terms
          end do
        end if
c end loop over Fourier components
      end do
c copy results back
      k = ( jrf - 1 ) * na + 1
      l = krf * na
      do i = k, l
        grdfld( i, itype ) = w( i )
      end do
c return local workspace
      deallocate ( w )
      deallocate ( gr )
      return
      end
