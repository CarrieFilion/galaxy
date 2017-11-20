      real*8 function ajrtab( E, L )
c  Copyright (C) 2017, Jerry Sellwood
      use aarrays
      implicit none
c this routine generates and then uses a look up table
c   uses CMLIB routines DB2INK and DB2VAL
c
c calling arguments
      real*8 L, E
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      real*8 actj1, db2val, Lzmax, splint2
c
c local arrays
      real*8, allocatable, save :: ajr(:,:)
      real*8, allocatable, save :: Bc1(:)
      real*8, allocatable, save :: Bc2(:)
      real*8, allocatable, save :: Eval(:)
      real*8, allocatable, save :: lamda(:)
      real*8, allocatable, save :: Lmx(:)
      real*8, allocatable, save :: Lvlu(:)
      real*8, allocatable, save :: tE(:)
      real*8, allocatable, save :: tL(:)
      real*8, allocatable, save :: w(:)
c
c local variables
      integer iE, ifail, iL, iuse, jcomp, lw, n, nE, nL
      real rm
      real*8 ajmax, E1, Emxel, Emnel, fp, L1, Lm
      save Emxel, Emnel, iuse, lw
c      parameter ( nE = 101, nL = 51, n = nE + nL )
      parameter ( nE = 201, nL = 101, n = nE + nL )
      include 'inc/pi.f'
c
      data iuse / 0 /
c
      if( iusead .eq. 0 )then
c allocate space on first call only
        if( iuse .eq. 0 )then
          allocate ( ajr( nL, nE ) )
          allocate ( Bc1( nE + 4 ) )
          allocate ( Bc2( nE * nL ) )
          allocate ( Eval( nE ) )
          allocate ( lamda( nE + 4 ) )
          allocate ( Lmx( nE ) )
          allocate ( Lvlu( nL ) )
          allocate ( tE( nE + 4 ) )
          allocate ( tL( nL + 4 ) )
          lw = nE * nL + 8 * ( max( nE, nL ) + 1 )
          allocate ( w( lw ) )
          iuse = 1
        end if
        print *, 'AJRTAB: building a table for actj1'
c use the largest radial range that will be needed
        jcomp = icomp
        rm = rmax
        call cutoff( sngl( arad( nradad ) ) )
        Emxel = Emaxe
        Emnel = Emine
        ajmax = 0
c choose evals - nominally y
        do iE = 1, nE
          E1 = dble( iE - 1 ) / dble( nE - 1 )
          E1 = 1. - cos( .5 * pi * E1 )
          Eval( iE ) = E1
          E1 = Emnel + E1 * ( Emxel - Emnel )
          Lm = Lzmax( E1 )
          Lmx( iE ) = Lm
c choose Lvalues - nominally x
          do iL = 1, nL
            if( iE .eq. 1 )then
              L1 = dble( iL - 1 ) / dble( nL - 1 )
              L1 = 1. - cos( .5 * pi * L1 )
              Lvlu( iL ) = L1
            else
              L1 = Lvlu( iL )
            end if
            L1 = L1 * Lm
            ajr( iL, iE ) = actj1( E1, L1 )
            ajmax = max( ajmax, ajr( iL, iE ) )
          end do
        end do
c fit spline interpolant to 2-D surface
        ifail = 0
        call db2ink( Lvlu, nL, Eval, nE, ajr, nL, 4, 4, tL, tE, Bc2,
     +               w, ifail )
        if( ifail .ne. 1 )then
          print *, 'IFLAG =', ifail
          call crash( 'AJRTAB', 'DB2INK failed' )
        end if
        iusead = 1
c fit spline to Lmax( E ) table
        E1 = Eval( 1 )
        Lm = splint2( Eval, Lmx, nE, E1, lamda, Bc1, .true. )
        print *, 'AJRTAB: table now ready'
c reset energy bounds
        call cutoff( rm )
        icomp = jcomp
        icmp = iccmp( icomp )
      end if
c
c interpolate Lm from table
      E1 = ( E - Emnel ) / ( Emxel - Emnel )
      E1 = max( E1, 0.d0 )
      E1 = min( E1, 1.d0 )
c ensure points on the boundary are deemed to be inside
      fp = E1 - Eval( nE )
      if( ( fp .ge. 0.d0 ) .and. ( fp .lt. 1.d-6 ) )E1 = Eval( nE )
      if( fp .ge. 1.d-6 )then
        print *, 'E1 out of range in ajrtab', E1, Eval( nE )
        print *, Emnel, E, Emxel
      end if
      Lm = splint2( Eval, Lmx, nE, E1, lamda, Bc1, .false. )
      L1 = abs( L )
      if( Lm .gt. 0.d0 )L1 = L1 / Lm
c      fp = L1 - Lvlu( nL )
c      if( ( fp .ge. 0.d0 ) .and. ( fp .lt. 1.d-6 ) )L1 = Lvlu( nL )
      E1 = max( E1, 0.d0 )
      E1 = min( E1, Eval( nE ) )
      L1 = max( L1, 0.d0 )
      L1 = min( L1, Lvlu( nL ) )
c look up ajr value in 2-D surface fit
      ajrtab = db2val( L1, E1, 0, 0, tL, tE, nL, nE, 4, 4, Bc2, w )
      if( ajrtab .eq. 0.d0 )call crash( 'AJRTAB', 'DB2VAL failed' )
      ajrtab = max( ajrtab, 0.d0 )
      return
      end
