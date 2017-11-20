      program tipsy
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c to convert a TIPSY file to a .pcs file
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      integer copenr, cread, lnblnk
      logical gtlogl
c
c local array
      integer nw
      real, allocatable :: w( : )
      real, allocatable :: gas( :,: )
      real, allocatable :: dark( :,: )
      real, allocatable :: star( :,: )
      real, allocatable :: outb( : )
      integer inp( 8 )
c
c local variables
      character fname*50, a*1, aw*4
      integer handle, tbytes
      integer i, it, j, k, l, n, ndim, ngas, ndark, nstar
      logical lswap
      real lunit, munit, pm, rt, rw, t( 2 ), vunit, xmn, xmx, ymn, ymx
      real*8 mdark, mgas, mstar, time
      equivalence ( it, rt ), ( rw, aw ), ( t, time )
c
    1 print *, 'Mass and velocity scaling depends on box size'
      call gtintg( '25 or 50 Mpc box?', i )
      if( i .eq . 25 )then
c unit conversions to solar masses and km/s for cosmo25cmb files
        munit = 2.310e15
        vunit = 630.28
      else if( i .eq. 50 )then
c unit conversions to solar masses and km/s for cosmo50cmb files
        munit = 1.84793e16
        vunit = 1261.05
      else
        go to 1
      end if
c set lunit in kpc
      lunit = i * 1000
c
c open file
      read '( a )', fname
      handle = copenr( fname( 1:i ) )
      lswap = gtlogl( 'Does data have the opposite endian rule' )
c read header
      j = 8
      tbytes = j * 4
      i = cread( handle, 4 * j, inp )
      if( lswap )then
        do i = 1, j
          call swap( inp( i ) )
        end do
      end if
      it = inp( 2 )
      t( 1 ) = rt
      it = inp( 1 )
      t( 2 ) = rt
      print *, 'time =', time
      n = inp( 3 )
      ndim = inp( 4 )
      ngas = inp( 5 )
      ndark = inp( 6 )
      nstar = inp( 7 )
      print *, '       n     ndim     ngas    ndark    nstar'
      print '( 6i9 )', n, ndim, ngas, ndark, nstar, inp( 8 )
      print *, 'If these values are nonsense then you probably' //
     +         ' guessed the wrong endian rule'
c
      nw = max( ngas * 12, ndark * 9, nstar * 11 )
      allocate ( w( nw ) )
c choose scaling to GALAXY units
      call set_scale
      lunit = lunit / unit_L
      munit = munit * 1.e-10 / unit_M
      vunit = vunit / unit_V
c read gas particles: m, x(3), v(3), rho, Y, H, Z, Phi
      it = ngas * 12
      tbytes = tbytes + it * 4
      if( it .gt. nw )call crash( 'MAIN', 'Input buffer too small 1' )
      i = cread( handle, 4 * it, w )
      if( lswap )then
        do i = 1, it
          call swap( w( i  ) )
        end do
      end if
      mgas = 0
      allocate ( gas( 7, ngas ) )
      k = 0
      do i = 1, ngas
        j = k + 1
        k = k + 12
        do l = 1, 3
          gas( l, i ) = lunit * w( l + j )
          gas( l + 3, i ) = vunit * w( l + 3 + j )
        end do
        gas( 7, i ) = munit * w( j )
        mgas = mgas + munit * w( j )
      end do
c read dark particles: m, x(3), v(3), eps, Phi
      it = ndark * 9
      tbytes = tbytes + it * 4
      if( it .gt. nw )call crash( 'MAIN', 'Input buffer too small 2' )
      i = cread( handle, 4 * it, w )
      if( lswap )then
        do i = 1, it
          call swap( w( i  ) )
        end do
      end if
      allocate ( dark( 7, ndark ) )
      k = 0
      rt = 0
      mdark = 0
      do i = 1, ndark
        j = k + 1
        k = k + 9
        do l = 1, 3
          dark( l, i ) = lunit * w( l + j )
          dark( l + 3, i ) = vunit * w( l + 3 + j )
        end do
        dark( 7, i ) = munit * w( j )
        mdark = mdark + munit * w( j )
      end do
c read star particles: m, x(3), v(3), Z, tform, eps, Phi
      it = nstar * 11
      tbytes = tbytes + it * 4
      if( it .gt. nw )call crash( 'MAIN', 'Input buffer too small 3' )
      i = cread( handle, 4 * it, w )
      print *, 'total number of bytes read', tbytes
      if( lswap )then
        do i = 1, it
          call swap( w( i  ) )
        end do
      end if
      allocate ( star( 7, nstar ) )
      k = 0
      mstar = 0
      do i = 1, nstar
        j = k + 1
        k = k + 11
        do l = 1, 3
          star( l, i ) = lunit * w( l + j )
          star( l + 3, i ) = vunit * w( l + 3 + j )
        end do
        star( 7, i ) = munit * w( j )
        mstar = mstar + munit * w( j )
      end do
      deallocate ( w )
      print *, 'mgas, mstar, mdark', sngl( mgas ), sngl( mstar ),
     +      sngl( mdark ), ' in units of 10^10 Msun'
c
c create .pcs file
      call new_unit( ndistf )
      ts = 1
      istep = 0
      call opnfil( ndistf, 'pcs', 'unformatted', 'new', 'seq', i )
      if( i .ne. 0 )then
        print *,'A file with this time stamp already exists'
        if( gtlogl( 'Do you want to overwrite it' ) ) then
          call opnfil( ndistf, 'pcs', 'unformatted', 'old', 'seq', i )
        else
          call crash( 'UNLOAD', 'Declined to overwite pcs file' )
        end if
      end if
c write header
      t( 1 ) = istep * ts
      print *, 'Creating .pcs file at time', t( 1 )
      k = 7
      l = 5000
c mean mass of a star particle
      pm = 1
      pertbn = .false.
      write( ndistf )nstar, ndark, ngas, k, l, t( 1 ), pm, pertbn
      allocate ( outb( k * l ) )
      j = 0
c output star particles
      if( nstar .gt. 0 )then
        do n = 1, nstar
          do i = 1, k
            outb( i + j ) = star( i, n )
          end do
          j = j + k
          if( j .ge. k * l )then
            write( ndistf )outb
            j = 0
          end if
        end do
      end if
c output dark matter particles
      if( ndark .gt. 0 )then
        do n = 1, ndark
          do i = 1, k
            outb( i + j ) = dark( i, n )
          end do
          j = j + k
          if( j .ge. k * l )then
            write( ndistf )outb
            j = 0
          end if
        end do
      end if
c output gas particles
      if( ngas .gt. 0 )then
        do n = 1, ngas
          do i = 1, k
            outb( i + j ) = gas( i, n )
          end do
          j = j + k
          if( j .ge. k * l )then
            write( ndistf )outb
            j = 0
          end if
        end do
      end if
      if( j .gt. 0 )write( ndistf )( outb( i ), i = 1, j )
      print *, 'created pcs file'
      end
