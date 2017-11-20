      subroutine c2grnm( direct )
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to create Green function files for the 2-D Cartesian grid
c
c It creates tables of the functions for x-accelerations, y-accelerations
c   and potentials, finding the Fourier transform of each and storing it
c   in a file.
c It will also save the direct (untransformed) data in a separate file if
c   requested by the calling argument being set to .true.
c If the files are already present, the routine simply reads the headers
c   and verifies that the values are consistent with the current run.
c
c This version uses single precision FFTPAK routines
      implicit none
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
      real, allocatable :: grn(:)
c
      integer ltrig
      real, allocatable :: trig(:)
c
      integer lwork
      real, allocatable :: w(:)
c
c local variables
      integer i, ifail, itype, ix, iy, j, jx, jy, lgrn, m, nx, ny
      real d2, x, y
      include 'inc/pi.f'
c
      if( .not. c2d )call crash( 'C2GRNM', 'Wrong type of method' )
c open files
      if( ngt .lt. 0 )call new_unit( ngt )
      call opnfil( ngt, 'grt', 'unformatted', 'old', 'seq', i )
      if( i .eq. 0 )then
c verify existing transformed file
        call grnhed( .true., .false., ifail )
        if( ifail .ne. 0 )call crash( 'C2GRNM',
     +                            'Wrong transformed Green fun file' )
        if( direct )then
          if( ngd .lt. 0 )call new_unit( ngd )
          call opnfil( ngd, 'grd', 'unformatted', 'old', 'seq', i )
          if( i .eq. 1 )call crash( 'C2GRNM',
     +                              'Direct Green fn file not found' )
c verify existing direct file
          call grnhed( .true., .true., ifail )
          if( ifail .ne. 0 )call crash( 'C2GRNM',
     +                                 'Wrong direct Green fun file' )
        end if
        print *, ' Old files found and checked in C2GRNM'
      else if( master )then
c should not create this file on the slaves as well as on the master
        print *, 'C2GRNM: Creating new Green fn file'
c allocate space
        lgrn = ( ngx + 1 ) * ( ngy + 1 )
        allocate ( grn( lgrn ) )
        lwork = 2 * lgrn
        allocate ( w( lwork ) )
        ltrig = 3 * max( ngx, ngy ) + 15
        allocate ( trig( ltrig ) )
c open file(s)
        call opnfil( ngt, 'grt', 'unformatted', 'new', 'seq', i )
        call grnhed( .false., .false., ifail )
        if( direct )then
          if( ngd .lt. 0 )call new_unit( ngd )
          call opnfil( ngd, 'grd', 'unformatted', 'new', 'seq', i )
          call grnhed( .false., .true., ifail )
        end if
c form x-forces, y-forces and then potentials
        do itype = 1, 3
          jx = 0
          if( itype .eq. 1 )jx = 1
          jy = 0
          if( itype .eq. 2 )jy = 1
c work over (x,y) field points
          m = 0
          do iy = jy, ngy
            y = dh( 2 ) * real( iy )
            do ix = jx, ngx
              x = dh( 1 ) * real( ix )
              m = m + 1
c Plummer softening is hard-wired!
              d2 = x * x + y * y + softl2
              if( itype .eq. 1 )grn( m ) = -x / d2**1.5
              if( itype .eq. 2 )grn( m ) = -y / d2**1.5
              if( itype .eq. 3 )grn( m ) = -1. / sqrt( d2 )
            end do
          end do
          if( direct )then
            write( ngd )( grn( i ), i = 1, m )
          end if
          nx = ngx + 1 - jx
          ny = ngy + 1 - jy
c Fourier transformation in x
          if( itype .eq. 1 )then
c Fourier sine transform of anti-symmetric function
            call rsinti( ngx - 1, trig )
            j = 1
            do iy = 1, ny
              call rsint( ngx - 1, grn( j ), trig )
              j = j + nx
            end do
          else
c Fourier cosine transform of symmetric function
            call rcosti( ngx + 1, trig )
            j = 1
            do iy = 1, ny
              call rcost( ngx + 1, grn( j ), trig )
              j = j + nx
            end do
          end if
c re-order array for Fourier transformation in y
          m = 0
          do ix = 1, nx
            j = ix
            do iy = 1, ny
              m = m + 1
              w( m ) = grn( j )
              j = j + nx
            end do
          end do
c Fourier transformation in y
          if( itype .eq. 2 )then
c find Fourier sine transform of anti-symmetric function
            call rsinti( ngy - 1, trig )
            j = 1
            do iy = 1, nx
              call rsint( ngy - 1, w( j ), trig )
              j = j + ny
            end do
          else
c find Fourier cosine transform of symmetric function
            call rcosti( ngy + 1, trig )
            j = 1
            do iy = 1, nx
              call rcost( ngy + 1, w( j ), trig )
              j = j + ny
            end do
          end if
c copy back
          m = 0
          do iy = 1, ny
            j = iy
            do ix = 1, nx
              m = m + 1
              grn( m ) = w( j )
              j = j + ny
            end do
          end do
c write out transformed version
          m = 0
          do iy = 1, ny - jy
            j = m + 1
            m = m + nx
            write( ngt )( grn( i ), i = j, m - jx )
          end do
        end do
        deallocate ( w )
        deallocate ( trig )
        deallocate ( grn )
c close and re-open file(s) to ensure that writes are completed
        close( ngt )
        call opnfil( ngt, 'grt', 'unformatted', 'old', 'seq', i )
        if( direct )then
          close( ngd )
          call opnfil( ngd, 'grd', 'unformatted', 'old', 'seq', i )
        end if
        write( no, * )' Green function creation complete'
        print *, 'C2GRNM: Green fn file ready'
      end if
      return
      end
