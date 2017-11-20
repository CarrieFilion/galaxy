      subroutine setgrd( pr )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up the grid to be used for the field determination
c   The input parameters will have been already read in
c Checks the grid size and sets the grid boundaries, interpolation scheme, etc
c
c Called from either GRDSET, CMPRUN or HEDREC
c
c calling argument
      logical pr
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      real grofu
      real*8 Gammaf, Gam1o1
c
c local variables
      character*8 str
      integer i, ifail, j, l, m, n
      logical lout, firstc
      real a, amax
      real*8 a1, a2, a3, deltam0, knl
      save firstc
      include 'inc/pi.f'
      data firstc / .true. /
c
      lout = master .and. pr
c set default value for all methods
      rinh( jgrid ) = 0
c estimate outer boundary
      amax = 1.1 * rmax
      if( amax .eq. 0. )amax = 10
      if( noslfg )then
c define some phony grid boundaries
        rgrid( jgrid ) = amax
        zm( jgrid ) = rgrid( jgrid )
c
      else if( p2d .or. p3d )then
c
c ensure that grid contains an even number of points in angle
        na = 2 * ( na / 2 )
c ensure that grid contains an suitable number of vertical mesh points
        if( p3d )then
          i = 0
          do while ( mod( ngz, 3 ) .eq. 0 )
            i = i + 1
            ngz = ngz / 3
          end do
          j = 0
          do while ( mod( ngz, 5 ) .eq. 0 )
            j = j + 1
            ngz = ngz / 5
          end do
          if( ngz .ne. 1 )then
            if( master )write( no, '( ''ngz ='', i4, '' x 3**'',' //
     +            ' i1, '' x 5**'', i1 )' )ngz, i, j
            call crash( 'SETGRD', 'ngz must be a multiple of 3 & 5' )
          end if
          ngz = 3**i * 5**j
        else
          ngz = 1
        end if
c grid size and Fourier harmonics
        if( lout )then
          if( p2d )then
            write( no, '( //'' Mesh size'', 2i6 )' )nr( jgrid ), na
          else
            write( no, '( //'' Mesh size'', 3i6 )' )nr( jgrid ), na,
     +                                              ngz
          end if
        end if
        mesh( jgrid ) = nr( jgrid ) * na * ngz
        mmax = na / 2 + 1
        alpha = 2. * pi / real( na )
        ng = min( ng, mmax )
        ng = min( ng, mmax )
        if( lout )write( no,
     + '( '' Sectoral harmonics 0 -'', i3, ''  will be used'' )' )ng - 1
c grid spacing exponent
        gamma = min( gamma, 1. )
        logrid = abs( 1. - gamma ) .lt. 0.01
        if( logrid )then
          gamma = 1
          ibeta = 0
          beta = -1
          if( lout )write( no, * )
     +                           ' Spacing of grid rings is logarithmic'
        else
          if( lout )write( no, '( 5x, ''gamma ='', f6.3 )' )gamma
          beta = 1. / ( 1. - gamma )
          ibeta = 1. - gamma
          if( lout )write( no,
     +     '('' Power law spacing of grid rings with exponent'', f8.3 )'
     +                                                             )beta
        end if
c inner edge of grid
        if( stdpol )then
          hole = .false.
        else
          if( logrid )then
            hole = uoffs .le. 0.
          else
            uoffs = max( uoffs, 0. )
            hole = uoffs .gt. 0.
          end if
        end if
        rinh( jgrid ) = grofu( 0. )
c set default for wholeg
        wholeg = .false.
c linear interpolation only for p2d
        if( p2d )jmass = 2
        if( p3d )then
          if( jmass .eq. 2 )then
            if( lout )write( no, * )'Linear interpolation'
          else if( jmass .eq. 3 )then
            if( lout )write( no, * )'Quadratic interpolation'
          else
            call crash( 'SETGRD', 'Invalid interpolation scheme' )
          end if
        end if
c radial grid boundary
        rgrid( jgrid ) = max( rgrid( jgrid ),
     +      grofu( real( nr( jgrid ) ) - .5 * real( jmass ) ) )
        if( lout )write( no,
     +           '( '' Inner and outer edges of grid are at'', f8.2,' //
     +             ' '' and'', f10.4, '' grid units'' )' )rinh( jgrid ),
     +                                                    rgrid( jgrid )
        if( ( .not. stdpol ) .and. ( .not. hole ) )then
          if( lout )write( no, '( '' There are'', i3,' //
     +               '  '' grid rings between r=0 and r='', f7.3 )' )
     +               nint( uoffs ) - 1, grofu( uoffs )
        end if
c vertical grid boundary
        if( p3d )then
          a = .5 * dzg *
     +                         real( ngz + 1 - jmass )
          if( ncode .eq. igrid( 1 ) )then
            zm( jgrid ) = a
          else
            zm( jgrid ) = max( a, zm( jgrid ) )
          end if
          if( lout )write( no,
     +'( '' Maximum z distance of grid is'', f8.2, '' grid units'' )' )a
        end if
c softening length
        softl2 = softl * softl
        if( lout )then
          write( no,
     +            '( '' Softening length is'', f10.6, '' grid units'' )'
     +                                                   )softl
          if( tsoft .eq. 1 )then
            write( no, * )'Plummer softening rule used'
          else if( tsoft .eq. 2 )then
            write( no, * )'Cubic spline softening rule used'
          end if
        end if
        if( ( tsoft .lt. 1 ) .or.
     +      ( tsoft .gt. 2 ) )then
          if( master )print *, 'tsoft =', tsoft
          call crash( 'SETGRD', 'Unrecognized softening rule' )
        end if
c set constants to be used in interpolation
        thmax = na
c allow for subdivision into sectors
        if( nsect .le. 0 )call crash( 'SETGRD', 'nsect not set' )
        thmax = thmax / real( nsect )
c check for sensible grid cell sizes in relation to softl
        a = nr( jgrid ) - 2
        a = rgrid( jgrid ) - grofu( a )
        if( a .gt. 3. * softl )then
          if( master )then
            print *
            write( str, '( f8.4 )' )a / softl
            print *, 'Warning: outermost cell sizes are' // str //
     +                         ' times larger than the softening length'
            print *, 'When this ratio >> 3, results may be unduly '
     +                                       // 'influenced by the grid'
            print *, 'Consider:'
            print *, '        increasing the softening length,'
            print *, ' and/or increasing the number of spokes,'
            print *, ' and/or reducing the number of rings'
          end if
          if( a / softl .gt. 5. )call crash( 'SETGRD',
     +                          'Unsuitable choice of grid parameters' )
          if( master )print *
        end if
c
      else if( c3d )then
c
c grid size
        if( lout )write( no,
     +               '( /'' Mesh size'', 3i6 )' )ngx, ngy, ngz
        ngxy = ngx * ngy
        mesh( jgrid ) = ngx * ngy * ngz
c linear interpolation is hard wired
        jmass = 2
        nplanes = ngz - jmass - 1
c coordinate boundaries in mesh spaces
        xm = .5 * real( ngx - jmass - 1 )
        ym = .5 * real( ngy - jmass - 1 )
        if( ncode .eq. igrid( 1 ) )then
          zm( jgrid ) = .5 * real( ngz - jmass - 1 )
        else
          a = .5 * real( ngz - jmass - 1 )
          zm( jgrid ) = max( a, zm( jgrid ) )
        end if
        a = min( xm, ym )
        rgrid( jgrid ) = max( rgrid( jgrid ), a )
c
      else if( c2d )then
c
c grid size
        if( lout )write( no, '( //'' Mesh size'', 3i6 )' )ngx, ngy
        ngxy = ngx * ngy
        mesh( jgrid ) = ngxy
c linear interpolation is hard wired
        jmass = 2
c assume square grid cells for now
        dh( 1 ) = 1
        dh( 2 ) = 1
c coordinate boundaries in mesh spaces
        xm = .5 * real( ngx - jmass - 1 )
        ym = .5 * real( ngy - jmass - 1 )
        a = min( xm, ym )
        rgrid( jgrid ) = max( rgrid( jgrid ), a )
c softening length
        softl2 = softl * softl
        if( lout )write( no, '( '' Softening length is'', f10.6,' //
     +                             ' '' mesh spaces'' )' )softl
c
      else if( sf2d .or. sf3d )then
c
        if( basset .eq. 'ablj' )then
          if( lout )write( no,
     +               '( ''Abel-Jacobi set'', i3, '' selected'' )' )basis
          softl = 0
          if( sf3d )call crash( 'SETGRD', 'Impossible basis' )
        else if( basset .eq. 'bess' )then
          if( lout )write( no, '( ''Bessel functions with deltak ='','
     +                              //  ' f6.3, '' selected'' )' )deltak
          if( lout )write( no, '( '' Softening length is'', f10.6,' //
     +                              ' '' grid units'' )' )softl
        else if( basset .eq. 'lgsp' )then
          if( lout )write( no, '( ''Log spiral functions with'' ' //
     +                   ' '' dalpha ='', f6.3, '' selected'' )' )deltak
          softl = 0
        else
          call crash( 'SETGRD', 'Unrecognized basis' )
        end if
        if( nsect .le. 0 )call crash( 'SETGRD', 'nsect not set' )
c determine min & max of n & m - needed when called from HEDREC
c    and doesn't hurt when called from GRDSET
        minm = 100
        maxm = 0
        minn = 100
        maxn = -100
        do i = 1, lastf
          minm = min( minm, msel( i ) )
          minn = min( minn, nsel( i ) )
          maxm = max( maxm, msel( i ) )
          maxn = max( maxn, nsel( i ) )
        end do
c set Fourier filter - not used except to report in hedrec
        ng = maxm + 1
        do i = 1, ng
          lg( i ) = .true.
        end do
        do i = 1, lastf
          j = msel( i ) + 1
          lg( j ) = .false.
        end do
        if( lout )then
          write( no,
     +              '( '' Number of basis functions used'', i4 )' )lastf
          write( no, * )'Radial (n) and azimuthal (m) basis functions'
          write( no, '( 3x, a )' )'minm maxm minn maxn'
          write( no, '( 2x, 4( 1x, i4 ) )' )minm, maxm, minn, maxn
        end if
        if( ( basset .eq. 'lgsp' ) .and.
     +      ( minm .eq. 0 ) )call crash( 'SETGRD',
     +                       'Axisymmetric log spirals not programmed' )
c number of planes
        if( sf3d )then
          ngz = max( ngz, 1 )
          a = dzg * real( ngz / 2 )
          if( lout )then
            write( no, '( 3x, a, i4 )')'number of planes ', ngz
            write( no, '( '' Maximum z distance of grid is'', f8.2' //
     +                                          ', '' grid units'' )' )a
          end if
c allow for predefined zm
          if( ncode .eq. igrid( 1 ) )then
            zm( jgrid ) = a
          else
            zm( jgrid ) = max( a, zm( jgrid ) )
          end if
        else
c one plane
          ngz = 1
        end if
        softl2 = softl * softl
c if no axisymmetric basis functions were chosen then an analytic
c   axisymmetric force is required
        fixrad = .true.
        do i = 1, lastf
          fixrad = fixrad .and. ( msel( i ) .ne. 0 )
        end do
        if( fixrad )then
          if( lout )write( no, * )' Using an analytic function' //
     +                    ' for the axisymmetric component of the force'
        end if
c set radius parameters for consistency with other methods
        rgrid( jgrid ) = max( maxmaxr, rgrid( jgrid ) )
        nr( jgrid ) = max( nr( jgrid ), 1 )
c number of coefficients
        ngx = 2 * maxm + 1
        ngy = maxn + 1
        ngxy = ngx * ngy
        mesh( jgrid ) = ngxy * ngz
c linear interpolation between planes
        jmass = 2
c
      else if( p3a )then
c
c ensure that grid contains an odd number of points in z
        ngz = 2 * ( ngz / 2 ) + 1
        if( lout )write( no,
     +            '( //'' Mesh size'', 2i6 )' )nr( jgrid ), ngz
        mesh( jgrid ) = nr( jgrid ) * ngz
        nsect = 1
c radial grid spacing exponent
        gamma = min( gamma, 1. )
        if( lout )write( no, '( 5x, ''gamma ='', f6.3 )' )gamma
        beta = 1. / ( 1. - gamma )
        ibeta = 1. - gamma
        if( lout )write( no, '( '' Power law spacing of grid rings'' '
     +                           //  ' '' with exponent'', f8.3 )' )beta
c grid boundaries
        rgrid( jgrid ) = max( rgrid( jgrid ),
     +                        grofu( real( nr( jgrid ) - 1 ) ) )
        a = dzg * real( ngz / 2 )
        if( ncode .eq. igrid( 1 ) )then
          zm( jgrid ) = a
        else
          zm( jgrid ) = max( a, zm( jgrid ) )
        end if
        if( lout )then
          write( no, '( '' Inner and outer edges of radial'' ' //
     +               ' '' grid are at'', f8.2, '' and'', f10.4, ' //
     +               ' '' grid units'' )' )grofu( 0. ), rgrid( jgrid )
          write( no, '( '' Maximum z distance of grid is'', f8.2,' //
     +                                 ' '' grid units'' )' )zm( jgrid )
        end if
c force type instruction
        if( lout )then
          if( fixrad )then
             write( no, * )'Radial forces held fixed'
          else
             write( no, * )'Radial forces active'
          end if
        end if
c linear interpolation between planes
        jmass = 2
c
      else if( s3d )then
c
        if( firstc )then
c spread points logarithmically in radius
          if( s3rad( nr( jgrid ) ) .le. 0. )call crash( 'SETGRD',
     +                                     'grid outer radius not set' )
          s3dex = log10( s3rad( nr( jgrid ) ) + 1. ) /
     +                                        real( nr( jgrid ) - 1 )
          if( nr( jgrid ) .gt. mr1 )then
            if( lout )print *, nr( jgrid ), mr1
            call crash( 'SETGRD', 'Radial arrays too small' )
          end if
          do i = 1, nr( jgrid )
            s3rad( i ) = grofu( real( i - 1 ) )
          end do
          firstc = .false.
        else
c check that this is no larger than the first
          do m = 1, jgrid - 1
            if( igrid( m ) .eq. 6 )then
              if( nr( jgrid ) .gt. nr( m )
     +     )call crash( 'SETGRD', 'first s3d grid must be the largest' )
            end if
          end do
        end if
c grid boundaries
        rgrid( jgrid ) = max( rgrid( jgrid ), s3rad( nr( jgrid ) ) )
        if( ncode .eq. igrid( 1 ) )then
          zm( jgrid ) = s3rad( nr( jgrid ) )
        else
          zm( jgrid ) = max( s3rad( nr( jgrid ) ), zm( jgrid ) )
        end if
        if( lout )then
          write( no, '( / i6, '' radial grid points.  Outermost'' '
     +       // ' '' is at'', f8.2 )' )nr( jgrid ), s3rad( nr( jgrid ) )
          write( no, '( ''Force terms up to l='', i2 )' )s3lmax
        end if
c linear interpolation between shells
        jmass = 2
c check space
        s3ntm = ( ( s3lmax + 1 ) * ( s3lmax + 2 ) ) / 2
        if( s3ntm .gt. s3mtm )call crash( 'SETGRD',
     +                        'parameter s3mtm in / grids / too small' )
        mesh( jgrid ) = 4 * s3ntm * nr( jgrid )
c
      else if( scf )then
c
        jmass = 1
c check space
        if( ( nbas .gt. nbasmx ) .or. ( s3lmax .gt. lbasmx ) )then
          if( lout )then
            print *, 'max allowed nbas, s3lmax', nbasmx, lbasmx
            print *, '  requested nbas, s3lmax', nbas, s3lmax
          end if
          call crash( 'SETGRD', 'arrays too small in / scfvars /' )
        end if
        s3ntm = ( ( s3lmax + 1 ) * ( s3lmax + 2 ) ) / 2
        if( s3ntm .gt. s3mtm )call crash( 'SETGRD',
     +                        'parameter s3mtm in / grids / too small' )
        mesh( jgrid ) = 2 * s3ntm * ( nbas + 1 )
c boundaries of volume for plotting and analysis
        if( rgrid( jgrid ) .le. 0. )then
          if( lout )print *, 'Nonsense rgrid', rgrid( jgrid )
          call crash( 'SETGRD', 'Nonsense rgrid' )
        end if
        if( ncode .eq. igrid( 1 ) )then
          zm( jgrid ) = rgrid( jgrid )
        else
          zm( jgrid ) = max( rgrid( jgrid ), zm( jgrid ) )
        end if
        if( lout )write( no,
     + '( ''Force terms up to n='', i2, '' and l='', i2 )' )nbas, s3lmax
        ifail = 0
        if( basset .eq. 'hern' )then
c initialize expansion coefficients - Hernquist basis
          do n = 0, nbas
            do l = 0, s3lmax
              knl = 0.5 * n * ( n + 4 * l + 3 ) +
     +                        ( l + 1 ) * ( 2 * l + 1 )
              a1 = n + 1
              a2 = 2. * l + 1.5
              a3 = n + 4 * l + 3
              anltilde( n, l ) = -2**( 8 * l + 6 ) *
     +                            ( n + 2 * l + 1.5 ) * Gam1o1( a1, a3 )
              anltilde( n, l ) = anltilde( n, l ) *
     +                        Gammaf( a2, ifail )**2 / ( 4. * pi * knl )
            end do
          end do
          do l = 0, s3lmax
            twoalpha( l ) = 2 * ( 2 * l + 1.5 )
          end do
        else if( basset .eq. 'plum' )then
c initialize expansion coefficients - Plummer (Clutton-Brock) basis
          do n = 0, nbas
            do l = 0, s3lmax
              knl = 4 * n * ( n + 2 * l + 2 ) +
     +                      ( 2 * l + 1 ) * ( 2 * l + 3 )
              a1 = n + 1
              a2 = l + 1
              a3 = n + 2 * l + 2
              anltilde( n, l ) = -2**( 4 * l + 4 ) *
     +                                  ( n + l + 1 ) * Gam1o1( a1, a3 )
              anltilde( n, l ) = anltilde( n, l ) *
     +                             Gammaf( a2, ifail )**2 / ( pi * knl )
            end do
          end do
          do l = 0, s3lmax
            twoalpha( l ) = 2 * ( l + 1 )
          end do
        else
          call crash( 'SETGRD', 'Unknown basis for SCF method' )
        end if
c remaining constants are the same for both bases
        do l = 0, s3lmax
          do m = 0, l
            deltam0 = 2.
            if( m .eq. 0 )deltam0 = 1.
            coeflm( l, m ) = ( 2. * l + 1. ) * deltam0 *
     +                  Gam1o1( dble( l - m + 1 ), dble( l + m + 1 ) )
          end do
        end do
        do n = 1, nbas
          c3( n ) = 1.0 / ( n + 1.0 )
          do l = 0, s3lmax
            c1( n, l ) = 2.0 * n + twoalpha( l )
            c2( n, l ) = n - 1.0 + twoalpha( l )
          end do
        end do
        lskip = 1
c        if( zerood .or. zeroev )lskip = 2
        lmin = 0
c        if( zeroev )lmin = 1
c
      else if( dr3d .or. bht )then
c
c softening length
        softl2 = softl * softl
        if( lout )write( no,
     +             '( '' Softening length is'', f10.6 )' )softl
c define some phony grid boundaries
        rgrid( jgrid ) = amax
        zm( jgrid ) = rgrid( jgrid )
c Barnes-Hut tree method has to be 3D
        if( bht .and. ( ndimen .ne. 3 ) )then
          print *, ndimen, 'D version of '
          call crash( 'SETGRD', 'Barnes-Hut tree scheme not available' )
        end if
c
      else
        call crash( 'SETGRD', 'Unrecognized method' )
      end if
c set flag that gravitational fields have not yet been calculated
      isfld = -1000
      return
      end
