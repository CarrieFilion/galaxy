      real*8 function satab( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c this routine generates and then uses a look up table
c   uses NAG routines E02CAF and E02CBF
c
c calling arguments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      include 'inc/orbval.f'
c
c externals
      real*8 db2val, Lzmax, Lzmin, sarea
c, rcirc, zmax
c
c local arrays
      integer nE, nLz
      parameter ( nE = 51, nLz = 11 )
      real*8, allocatable :: area(:,:), Bc2(:), Eval(:), Lzv1(:)
      real*8, allocatable :: tE(:), tL(:), work(:)
      save Bc2, tE, tL, work
c
c local variables
      integer iE, ifail, iLz, iuse, lwork
      real*8 E1, Elast, Lz1, Lzl, Lzm
c, r, aa
      save Elast, Lzl, Lzm, iuse
      include 'inc/pi.f'
c
      data iuse / 0 /, Elast / -100. /
c
      if( iuse .eq. 0 )then
        if( master )print *, 'SATAB: building a table'
        allocate ( Eval( nE ) )
        allocate ( Lzv1( nLz ) )
        allocate ( area( nLz, nE ) )
c choose Lzvals - concentrated towards the ends of the range
        do iLz = 1, nLz
          Lz1 = real( iLz - 1 ) / real( nLz - 1 )
          Lz1 = .5 * ( 1.d0 - cos( pi * Lz1 ) )
          Lz1 = max( Lz1, 0.d0 )
          Lz1 = min( Lz1, 1.d0 )
          Lzv1( iLz ) = Lz1
        end do
c choose Evals - concentrated towards the ends of the range
        do iE = 1, nE
          E1 = real( iE - 1 ) / real( nE - 1 )
          E1 = .5 * ( Emaxe + Emine -
     +                              ( Emaxe - Emine ) * cos( pi * E1 ) )
          E1 = max( E1, Emine )
          E1 = min( E1, Emaxe )
          Eval( iE ) = E1
          Lzm = Lzmax( E1 )
          Lzl = Lzmin( E1 )
          Lzm = max( Lzm, Lzl + 1.d-8 )
c work over Lz values
          do iLz = 1, nLz
            Lz1 = Lzv1( iLz )
            Lz1 = Lzm * Lz1 + Lzl * ( 1.d0 - Lz1 )
            Lz1 = max( Lz1, Lzl )
            Lz1 = min( Lz1, Lzm )
            area( iLz, iE ) = sarea( E1, Lz1 )
          end do
        end do
c fit a 2D cubic spline
        allocate ( tL( nLz + 4 ) )
        allocate ( tE( nE + 4 ) )
        allocate ( Bc2( nLz * nE ) )
        lwork = nLz * nE + 8 * ( max( nLz, nE ) + 1 )
        allocate ( work( lwork ) )
        ifail = 0
        call db2ink( Lzv1, nLz, Eval, nE, area, nLz, 4, 4, tL, tE,
     +               Bc2, work, ifail )
        if( ifail .gt. 1 )then
          if( master )print *, 'iflag =', ifail
          call crash( 'SATAB', 'DB2INK failed' )
        end if
        if( master )print *, 'SATAB: table now ready'
        iuse = 1
      end if
c general potential
      if( E .ne. Elast )then
        Lzl = Lzmin( E )
        Lzm = Lzmax( E )
        Elast = E
      end if
c particle at rest in the centre
      if( Lzm .eq. 0. )then
        satab = 0
      else
c allow for points outside fitted range only because of round-off
        E1 = min( E, Emaxe )
        E1 = max( E1, Emine )
        if( abs( E - E1 ) .gt. 1.d-6 )then
          print *, 'Attempt to call db2val with', sngl( Lz ), sngl( E )
          print *, 'bounds are', sngl( Lzl ), sngl( Lzm ),
     +                                      sngl( Emine ), sngl( Emaxe )
          call crash( 'SATAB', 'point outside fitted area' )
        end if
c interpolate value from table
        Lz1 = ( abs( Lz ) - Lzl ) / ( Lzm - Lzl )
        satab = db2val(
     +                 Lz1, E1, 0, 0, tL, tE, nLz, nE, 4, 4, Bc2, work )
      end if
      return
      end
