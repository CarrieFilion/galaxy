      real*8 function E0tab( a1, a2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c this routine generates and then uses a look up table
c   uses CMLIB routines DB2INK and DB2VAL
c
c calling arguments
      real*8 a1, a2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      real*8 aj1max, algrng2, db2val, Ej1j2, Emin, splint2
c
c local arrays
      integer na1, na2, n
      parameter ( na1 = 81, na2 = 81, n = na1 + na2 )
      real*8, allocatable, save :: a1mx(:,:)
      real*8, allocatable, save :: a1val(:,:)
      real*8, allocatable, save :: a2val(:,:)
      real*8, allocatable, save :: ac(:,:)
      real*8, allocatable, save :: alamda(:,:)
      real*8, allocatable, save :: Bc2(:,:)
      real*8, allocatable :: E0(:,:)
      real*8, allocatable, save :: Ec(:,:)
      real*8, allocatable, save :: Elamda(:,:)
      real*8, allocatable, save :: Emn(:,:)
      real*8, allocatable, save :: ta1(:,:)
      real*8, allocatable, save :: ta2(:,:)
      real*8, allocatable, save :: w(:)
c
c local variables
      character*6 fname
      integer i, ia1, ifail, ia2, iuse( 2 ), j, jcmp, juse, lw, unit
      logical lcmp, precalc, cuspy( 2 )
      real*8 aj1, aj1m, aj2, E1, E2, fp, Lzmx( 2 ), r, rm, rs
      save cuspy, iuse, juse, Lzmx
      include 'inc/pi.f'
c
      data iuse, juse / 0, 0, 0 /
c
      if( iuse( icomp ) .ne. icomp )then
        allocate ( E0( na1, na2 ) )
        if( juse .eq. 0 )then
          allocate ( a1mx( na2, ncomp ) )
          allocate ( a1val( na1, ncomp ) )
          allocate ( a2val( na2, ncomp ) )
          allocate ( ac( na2 + 4, ncomp ) )
          allocate ( alamda( na2 + 4, ncomp ) )
          allocate ( Bc2( na1 * na2, ncomp ) )
          allocate ( Ec( na2 + 4, ncomp ) )
          allocate ( Elamda( na2 + 4, ncomp ) )
          allocate ( Emn( na2, ncomp ) )
          allocate ( ta1( na1 + 4, ncomp ) )
          allocate ( ta2( na2 + 4, ncomp ) )
          lw = na1 * na2 + 8 * ( max( na1, na2 ) + 1 )
          allocate ( w( lw ) )
          juse = 1
        end if
c check component flags
        if( ( icomp .lt. 1 ) .or. ( icomp .gt. ncomp ) )then
          print *, 'icomp =', icomp
          call crash( 'E0TAB', 'Nonsense icomp' )
        end if
        i = 1
        do while ( i .ne. iccmp( icomp ) )
          i = i + 1
          if( i .gt. ncmp )call crash( 'E0TAB', 'Unrecognized iccmp' )
        end do
        if( i .ne. icmp )call crash( 'E0TAB', 'icmp/icomp mismatch' )
c set up single component if needed
        lcmp = cmprssd
        jcmp = icmp
        icmp = iccmp( icomp )
        rs = rmax
        call cutoff( rtrunc( icmp ) )
        j = ihalo0( icomp )
        icmp = jcmp
c recognize cuspy models
        cuspy( icomp ) = .false.
        cuspy( icomp ) = cuspy( icomp ) .or. ( ( j .eq. 9 ) .or.
     +             ( j .eq. 20 ) .or. ( j .eq. 25 ) .or. ( j .eq. 32 ) )
        if( ( ( j .ne. 4 ) .and. ( j .ne. 6 ) .and.
     +        ( j .ne. 16 ) ) .and. ( .not. cuspy( icomp ) )
     +                       )call crash( 'E0tab', 'Unrecognized halo' )
c attempt to open old .dat file
        fname = 'E0tab'
        j = 5
        if( ncomp .gt. 1 )then
          write( fname( 6:6 ), '( i1 )' )icomp
          j = 6
        end if
        call new_unit( unit )
        open( unit, file = fname( 1:j )//'.dat', form = 'unformatted',
     +        status = 'old', iostat = ifail )
        precalc = ifail .eq. 0
        r = rtrunc( icmp )
c check header
        if( precalc )then
          read( unit )rm, E1, E2
          fp = Emaxe - Emine
          precalc = ( abs( E1 - Emine ) / fp .lt. 1.d-6 ) .and.
     +              ( abs( E2 - Emaxe ) / fp .lt. 1.d-6 ) .and.
     +              ( abs( rm - r ) .lt. 1.d-6 )
          if( .not. precalc )then
            print *, E1, Emine, E1 - Emine
            print *, E2, Emaxe, E2 - Emaxe
            print *, rm, r, rm - r
          end if
        end if
        if( precalc )then
c checked out - read in precalculated data
          read( unit )( ( E0( i, j ), i = 1, na1 ), j = 1, na2 )
        else
c need to remake the file
          print *, 'New and old rmax', sngl( r ), sngl( rm )
          rm = r
          if( master )print
     +         '( ''E0TAB: Building a table for component'', i2 )', icmp
        end if
        close( unit )
        Lzmx( icomp ) = Lztmx( icmp )
c choose a2vals - nominally x
        do ia2 = 1, na2
          aj2 = dble( ia2 - 1 ) / dble( na2 - 1 )
          if( cuspy( icomp ) )then
            aj2 = aj2**2
            aj2 = 1. - cos( .5 * pi * aj2 )
          end if
          a2val( ia2, icomp ) = aj2
          aj2 = Lzmx( icomp ) * aj2
          aj1m = aj1max( aj2 )
          a1mx( ia2, icomp ) = aj1m
          Emn( ia2, icomp ) = Emin( aj2 )
c choose a1vals - nominally y
          do ia1 = 1, na1
            if( ia2 .eq. 1 )then
              aj1 = dble( ia1 - 1 ) / dble( na1 - 1 )
              if( cuspy( icomp ) )then
                aj1 = aj1**2
                aj1 = 1. - cos( .5 * pi * aj1 )
              end if
              a1val( ia1, icomp ) = aj1
            else
              aj1 = a1val( ia1, icomp )
            end if
            aj1 = aj1 * aj1m
            if( .not. precalc )E0( ia1, ia2 ) = ej1j2( aj1, aj2 )
          end do
        end do
        if( .not. precalc )then
          if( master )print *, 'E0TAB: Table now ready'
c save table values in a file
          fname = 'E0tab'
          j = 5
          if( ncomp .gt. 1 )then
            write( fname( 6:6 ), '( i1 )' )icomp
            j = 6
          end if
          open( unit, file = fname( 1:j )//'.dat', form = 'unformatted',
     +          status = 'unknown' )
          write( unit )rm, Emine, Emaxe
          write( unit )( ( E0( i, j ), i = 1, na1 ), j = 1, na2 )
          close( unit )
        end if
        ifail = 0
        call db2ink( a1val( 1, icomp ), na1, a2val( 1, icomp ), na2,
     +               E0, na1, 4, 4, ta1( 1, icomp ), ta2( 1, icomp ),
     +               Bc2( 1, icomp ), w, ifail )
        if( ifail .gt. 1 )then
          print *, 'iflag =', ifail
          call crash( 'E0TAB', 'DB2INK failed' )
        end if
        deallocate( E0 )
        if( cuspy( icomp ) )then
c fit spline interpolant to a1max curve
          aj2 = a2val( 1, icomp )
          aj1 = splint2( a2val( 1, icomp ), a1mx( 1, icomp ), na2, aj2,
     +       alamda( 1, icomp ), ac( 1, icomp ), .true. )
c fit spline interpolant to Emin curve
          aj2 = a2val( 1, icomp )
          aj1 = splint2( a2val( 1, icomp ), Emn( 1, icomp ), na2, aj2,
     +       Elamda( 1, icomp ), Ec( 1, icomp ), .true. )
        end if
c restore composite model
        if( lcmp )then
          cmprssd = .true.
          snglph = .false.
          call cutoff( sngl( rs ) )
        end if
        iuse( icomp ) = icomp
      end if
c interpolate value from table
      aj2 = abs( a2 ) / Lzmx( icomp )
      aj1 = max( a1, 0.d0 )
c ensure points on the boundary are deemed to be inside
      fp = aj2 - a2val( na2, icomp )
      if( ( fp .ge. 0.d0 ) .and.
     +    ( fp .lt. 1.d-6 ) )aj2 = a2val( na2, icomp )
c angular momentum cannot exceed the maximum in the uncompressed halo
      if( aj2 .gt. a2val( na2, icomp ) )then
        E0tab = Emaxo( icomp ) + 1
      else
        if( cuspy( icomp ) )then
          aj1m = splint2( a2val( 1, icomp ), a1mx( 1, icomp ), na2, aj2,
     +      alamda( 1, icomp ), ac( 1, icomp ), .false. )
        else
         aj1m = algrng2( a2val( 1, icomp ), a1mx( 1, icomp ), na2, aj2 )
        end if
        if( ( aj1m .gt. 1.d-6 ) .or. ( .not. cuspy( icomp ) ) )then
          if( aj1m .gt. 0.d0 )aj1 = aj1 / aj1m
          fp = aj1 - a1val( na1, icomp )
          if( ( fp .ge. 0.d0 ) .and.
     +         ( fp .lt. 1.d-3 ) )aj1 = a1val( na1, icomp )
c flag actions out of original range
          if( aj1 .lt. ta1( 1, icomp ) .or.
     +        aj1 .gt. ta1( na1 + 4, icomp ) .or.
     +        aj2 .lt. ta2( 1, icomp ) .or.
     +        aj2 .gt. ta2( na2 + 4, icomp ) )then
            E0tab = Emaxo( icomp ) + 1
          else
c interpolate value from spline fit to 2-D surface
            E0tab = db2val( aj1, aj2, 0, 0, ta1( 1, icomp ),
     +             ta2( 1, icomp ), na1, na2, 4, 4, Bc2( 1, icomp ), w )
          end if
c circular orbits
        else
          E0tab = splint2( a2val( 1, icomp ), Emn( 1, icomp ), na2, aj2,
     +      Elamda( 1, icomp ), Ec( 1, icomp ), .false. )
        end if
      end if
      return
      end
