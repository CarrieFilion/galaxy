      subroutine splfit
c  Copyright (C) 2015, Jerry Sellwood
c
c Creates and saves a cubic spline fit to azimuthally or spherically averaged
c   potential on the grid.  The spline parameters are calculated and saved in
c   natural units
c   uses NAG routines either E02BEF or E02DCF
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'inc/setup.f'
c
      integer iuseact
      real*8 El, Em, hlast, hm
      common / tabactn / hm, hlast, El, Em, iuseact
c
c local allocatable arrays
      real*8, allocatable :: Bc1(:), Bc2(:), pt(:), rt(:), sfw(:)
      real*8, allocatable :: tr(:), tz(:), zt(:)
c
c externals - azavpt works in grid units
      real azavpt
      real*8 gmassh
c
c externals
      logical gtlogl
      real grofu
      real*8 phihal, splint2
c
c local variables
      integer ip, ir, ifail, iz, j, jcmp, lr, lw, lz, n, uni
      logical firstc
      real dr, dz, rs, stype, zs
      real*8 fp, gm, s
      save firstc, uni
c
      data firstc / .true. /
c
      if( iusepf .ne. 0 )call crash( 'SPLFIT', 'Fit already made' )
      call switch( ngrid )
c allocate space
      lr = 0
      if( sphrod( icmp ) )then
        if( c3d )then
          lr = ngx / 2
          lz = ngz / 2
        else if( p3d .or. p3a )then
          lr = nr( jgrid )
          lz = ngz / 2
        else if( s3d )then
          lr = nr( jgrid )
          lz = lr
        end if
        allocate ( pt( lr * lz ) )
        allocate ( rt( lr ) )
        allocate ( zt( lz ) )
c set up input tables
        if( c3d )then
          n = 0
          do ir = 1, lr
            rs = ir - 1
            rt( ir ) = rs / lscale
            do iz = 1, lz
              zs = iz - 1
              if( ir .eq. 1 )zt( iz ) = zs / lscale
              n = n + 1
              pt( n ) = azavpt( rs, zs )
            end do
          end do
c mass distribution is axisymmetric - so p3d and p3a are equivalent
        else if( p3d .or. p3a )then
          n = 0
          do ir = 1, lr
            rs = grofu( real( ir - 1 ) )
            rt( ir ) = rs / lscale
            do iz = 1, lz
              zs = real( iz - 1 ) * dzg
              if( ir .eq. 1 )zt( iz ) = zs / lscale
              n = n + 1
              pt( n ) = azavpt( rs, zs )
            end do
          end do
        else if( s3d )then
          dr = rgrid( jgrid ) / real( lr - 1 )
          dz = zm( jgrid ) / real( lz - 1 )
          n = 0
          do ir = 1, lr
            rs = dr * real( ir - 1 )
            rt( ir ) = rs / lscale
            do iz = 1, lz
              zs = dz * real( iz - 1 )
              if( ir .eq. 1 )zt( iz ) = zs / lscale
              n = n + 1
              pt( n ) = azavpt( rs, zs )
            end do
          end do
        end if
      else
        if( s3d )then
          lr = nr( jgrid )
          allocate ( pt( lr ) )
          allocate ( rt( lr ) )
          dr = rgrid( jgrid ) / real( lr - 1 )
          lz = 1
          zs = 0
          do ir = 1, lr
            rs = dr * real( ir - 1 )
            rt( ir ) = rs / lscale
            pt( ir ) = azavpt( rs, zs )
          end do
        end if
      end if
      if( lr .eq. 0 )call crash( 'SPLFIT', 'Unrecognized grid' )
c add potential of any additional components
      if( ncmp .gt. 2 )then
        jcmp = icmp
        do icmp = 1, ncmp
          if( ( icmp .ne. jcmp ) .and. ( .not. disc( icmp ) ) )then
            n = 0
            do ir = 1, lr
              if( sphrod( jcmp ) )then
                do iz = 1, lz
                  s = sqrt( rt( ir )**2 + zt( iz )**2 )
                  n = n + 1
                  pt( n ) = pt( n ) + phihal( s )
                end do
              else
                s = rt( ir )
                n = n + 1
                pt( n ) = pt( n ) + phihal( s )
              end if
            end do
          end if
        end do
c adjust offset such that phi( rgrid, 0 ) = -GM/rgrid
        n = lr
        if( sphrod( jcmp ) )n = ( lr -1 ) * lz + 1
        gm = 0
        s = rgrid( jgrid ) / lscale
        do icmp = 1, ncmp
          if( .not. disc( icmp ) )then
            if( icmp .eq. jcmp )then
              gm = gm - pt( n ) * s
            else
              gm = gm + gmassh( s )
            end if
          end if
        end do
        fp = -gm / s - pt( n )
        n = 0
        do ir = 1, lr
          if( sphrod( jcmp ) )then
            do iz = 1, lz
              n = n + 1
              pt( n ) = fp + pt( n )
            end do
          else
            n = n + 1
            pt( n ) = fp + pt( n )
          end if
        end do
        icmp = jcmp
      end if
      if( sphrod( icmp ) )then
c fit a 2D cubic spline
        allocate ( tr( lr + 4 ) )
        allocate ( tz( lz + 4 ) )
        allocate ( Bc2( lr * lz ) )
        lw = lr * lz + 8 * ( max( lr, lz ) + 1 )
        allocate ( sfw( lw ) )
        ifail = 0
        call db2ink( zt, lz, rt, lr, pt, lz, 4, 4, tz, tr,
     +               Bc2, sfw, ifail )
        if( ifail .gt. 1 )then
          if( master )print *, 'iflag =', ifail
          call crash( 'SPLFIT', 'DB2INK failed' )
        end if
      else
c spherical mode - 1D spline
        allocate ( tr( lr + 4 ) )
        allocate ( Bc1( lr + 4 ) )
        s = .5 * rt( lr )
        fp = splint2( rt, pt, lr, s, tr, Bc1, .true. )
      end if
c save fit
      uni = -1
      if( firstc )then
        call opnfil( uni, 'pft', 'unformatted', 'new', 'seq', ifail )
        if( ifail .ne. 0 )then
          if( .not. gtlogl(
     +         '.pft file exists, do you want to over-write it?' ) )stop
          call opnfil( uni, 'pft', 'unformatted', 'old', 'seq', ifail )
        end if
        firstc = .false.
      else
        call opnfil( uni, 'pft', 'unformatted', 'old', 'seq', ifail )
      end if
      write( uni )lr + 4, lz + 4, ncmp, icmp, ncode
      do ip = 1, ncmp
        read( ctype( ip ), '( a4 )' )stype
        write( uni )imtyp( ip ), impar( ip ), idftyp( ip ),
     +              ( dfcns( j, ip ), j = 1, 3 ), cmpmas( ip ),
     +              fmass( ip ), rscale( ip ), rtrunc( ip ), disc( ip ),
     +              stype
      end do
      if( sphrod( icmp ) )then
        write( uni )( tr( ir ), ir = 1, lr + 4 ),
     +              ( tz( ir ), ir = 1, lz + 4 ),
     +              ( Bc2( ir ), ir = 1, lr * lz )
      else
         write( uni )( tr( ir ), ir = 1, lr + 4 ),
     +               ( Bc1( ir ), ir = 1, lr + 4 )
      end if
      close( uni )
c flag new table, also that action table needs to be recalculated
      initls = .true.
      iuseact = 0
      return
      end
