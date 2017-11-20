      real function zdisp( r, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Returns the vertical component of the velocity dispersion tensor of disc
c   particles at the point (r,z).  The value is that required to maintain
c   local vertical balance in the grid calculated vertical force field.  The
c   dispersion is determined by an integration of the local Jeans equation
c   [(4-36), p199 of Binney & Tremaine]
c The routine creates and stores a 2-D table of values on the first call
c   and subsequent calls return bi-linear interpolates from this table
c
c calling arguments and returned value are in model units, not grid units
      real r, z
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'inc/setup.f'
c
c externals
      real zthick
      real*8 gsigmt, rhorz
c
c local allocatable arrays
      real, allocatable :: fzmean(:,:), zset(:,:)
c
c local arrays
      integer n
      parameter ( n = 5 )
      real*8 absc( n ), wght( n )
c
c local variables
      integer i, ir, iz, iuse, jz, kz, npl, nz
      parameter ( nz = 51 )
      real dez, dz, fr, frcfac, fz, rs, zs
      real*8 a, b, rho, r2, z2
      include 'inc/pi.f'
      save dz, iuse, zset
c
      data iuse / 0 /
c
c skip if zero thickness disc
      zdisp = 0
      if( z0init( icmp ) .gt. 0. )then
c analytic Spitzer sheet
        if( iztyp( icmp ) .eq. 1 )then
          zdisp = sqrt( 2. * pi * gsigmt( dble( r ) ) * zthick( r ) )
c
        else
c
          call switch( 1 )
          npl = ngz - 1
          if( c3d )npl = nplanes
          if( s3d )then
            npl = 70
            dzg = .2 * zthick( 0. )
          end if
          if( npl .le. 1 )call crash( 'ZDISP', 'Nonsense npl value' )
c
          if( iuse .ne. icmp )then
            if( master )print *, 'ZDISP: Building a table for pop', icmp
c allocate space
            allocate ( fzmean( mradq + 1, npl ) )
            allocate ( zset( mradq + 1, nz ) )
c set radial spacing
            nradq = mradq
            r2 = max( rhole, dble( rtrunc( icmp ) ) / real( nradq ) )
            if( drq .le. 0.
     +        )drq = ( dble( rtrunc( icmp ) ) - r2 ) / real( nradq - 1 )
c tabulate mean vertical forces
            frcfac = 1. / ( lscale * ts**2 )
            do iz = 1, npl + 1
              z2 = dzg * ( iz - npl / 2 - 1 )
              zs = z2
              jz = abs( iz - npl / 2 - 1 ) + 1
              if( jz .gt. npl )call crash( 'ZDISP',
     +                                        'fzmean table too small' )
              if( jz .gt. 1 )then
                do ir = 1, nradq + 1
                  r2 = drq * real( ir - 1 ) * lscale
                  rs = r2
                  call mnfrz( rs, zs, fr, fz )
                  if( zs .lt. 0. )then
                    fzmean( ir, jz ) = fz
                  else
                    fzmean( ir, jz ) =
     +                           .5 * frcfac * ( fzmean( ir, jz ) - fz )
                  end if
                end do
              else
                do ir = 1, nradq + 1
                  fzmean( ir, jz ) = 0
                end do
              end if
            end do
c set up table spanning 7 scale heights
            dz = min( 7. * zthick( 0. ), zm( 1 ) )
            dz = dz / real( nz - 1 )
            do jz = 1, nz
              do ir = 1, nradq + 1
                r2 = drq * real( ir - 1 )
                z2 = dz * real( jz - 1 )
                rho = rhorz( r2, z2 )
                if( rho .gt. 0. )then
c piece-wise integration
                  zdisp = 0
                  kz = z2 * lscale / dzg + 1.
                  b = z2
                  do iz = kz, npl / 2
                    a = b
                    b = real( iz ) * dzg / lscale
                    z2 = a
                    if( rhorz( r2, z2 ) .gt. 0.d0 )then
c non-adaptive Gauss-Legendre quadrature
                      call GLqtab( a, b, n, wght, absc )
                      do i = 1, n
                        z2 = absc( i )
                        dez = z2 * lscale / dzg - real( iz - 1 )
                        zdisp = zdisp + wght( i ) * rhorz( r2, z2 ) *
     +                             ( fzmean( ir, iz ) * ( 1. - dez ) +
     +                                      fzmean( ir, iz + 1 ) * dez )
                      end do
                    end if
                  end do
                  zdisp = zdisp / rho
                  zset( ir, jz ) = sqrt( max( zdisp, 0. ) )
                else
                  zset( ir, jz ) = 0.
                end if
              end do
            end do
            iuse = icmp
            if( master )print *, 'ZDISP: Table ready'
c return local space
            deallocate ( fzmean )
          end if
c look up tabulated value
          rs = r / drq
          ir = rs + 1
          ir = min( ir, nradq )
          rs = rs - real( ir - 1 )
          zs = abs( z ) / dz
          jz = zs + 1.
          jz = min( jz, nz - 1 )
          zs = zs - real( jz - 1 )
          zdisp = zset( ir    , jz     ) * ( 1. - rs ) * ( 1. - zs ) +
     +            zset( ir + 1, jz     ) *     rs      * ( 1. - zs ) +
     +            zset( ir    , jz + 1 ) * ( 1. - rs ) *     zs +
     +            zset( ir + 1, jz + 1 ) *     rs      *     zs
        end if
      end if
      return
      end
