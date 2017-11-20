      real*8 function tautab( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c this routine generates and then uses a look up table
c    uses CMLIB routines DB2INK and DB2VAL
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
c externals
      real*8 akappa, algrng2, db2val, Lzmax, Lzmin, rcofe, tau, vcirc
c
c local arrays
      integer lwork, nE, nLz
      parameter ( nE = 51, nLz = 26 )
      real*8, allocatable :: Bc2(:), Eval(:), Lzv1(:)
      real*8, allocatable :: period(:,:), period1(:), tE(:), tL(:)
      real*8, allocatable :: work(:)
      save Bc2, Eval, period1, tE, tL, work
c
c local variables
      integer iE, ifail, iLz, iuse
      real*8 ak, E1, Elast, Lz1, Lzl, Lzm, r
      save ak, Elast, Lzl, Lzm, iuse
      include 'inc/pi.f'
c
      data iuse / 0 /, Elast / -100. /
c
      if( iuse .ne. icmp )then
c Mestel/Toomre/Zang disc
        if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'MTZ ' ) )then
          E1 = 1.
          Lzm = Lzmax( E1 )
          allocate ( Eval( nE ) )
          allocate ( period1( nE ) )
          do iLz = 1, nE
            Lz1 = Lzm * real( iLz - 1 ) / real( nE - 1 )
            Eval( iLz ) = Lz1
            period1( iLz ) = tau( E1, Lz1 )
          end do
        else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'ISOC' ) )then
c isochrone disc
          allocate ( Eval( nE ) )
          allocate ( period1( nE ) )
          do iE = 1, nE
            E1 = Emine + ( Emaxe - Emine ) * real( iE - 1 )
     +                                                  / real( nE - 1 )
            Eval( iE ) = E1
            E1 = max( E1, Emine )
            E1 = min( E1, Emaxe )
            r = rcofe( E1 )
            period1( iE ) = 2. * pi / akappa( r )
          end do
        else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'POWE' ) )then
c power law disc
          E1 = sign( 1., dfcns( 3, 1 ) )
          Lzm = Lzmax( E1 )
          allocate ( Eval( nE ) )
          allocate ( period1( nE ) )
          do iLz = 1, nE
            Lz1 = Lzm * real( iLz - 1 ) / real( nE - 1 )
            Eval( iLz ) = Lz1
            period1( iLz ) = tau( E1, Lz1 )
          end do
c Maclaurin/Freeman/Kalnajs disc is a const
        else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'MFK ' ) )then
          E1 = 1
c so is a uniform sphere
        else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'UNIS' ) )then
          E1 = 1
        else
c general potential
          if( master )print *, 'TAUTAB: building a table'
          allocate ( Eval( nE ) )
          allocate ( Lzv1( nLz ) )
          allocate ( period( nLz, nE ) )
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
            E1 = .5 * ( 1.d0 - cos( pi * E1 ) )
            E1 = E1 * Emaxe + ( 1.d0 - E1 ) * Emine
            E1 = max( E1, Emine )
            E1 = min( E1, Emaxe )
            Eval( iE ) = E1
            r = rcofe( E1 )
            Lzm = r * vcirc( r )
            ak = akappa( r )
            Lzl = Lzmin( E1 )
            if( Lzm - Lzl .lt. 1.d-8 )then
c since all orbits at this E1 must be closely circular, tau = 2pi/kappa
              do iLz = 1, nLz
                period( iLz, iE ) = 2. * pi
              end do
            else
c work over Lz values
              do iLz = 1, nLz
                Lz1 = Lzv1( iLz )
                Lz1 = Lzm * Lz1 + Lzl * ( 1.d0 - Lz1 )
                Lz1 = max( Lz1, Lzl )
                Lz1 = min( Lz1, Lzm )
                period( iLz, iE ) = tau( E1, Lz1 ) * ak
              end do
            end if
          end do
c fit a 2D cubic spline
          allocate ( tL( nLz + 4 ) )
          allocate ( tE( nE + 4 ) )
          allocate ( Bc2( nLz * nE ) )
          lwork = nLz * nE + 8 * ( max( nLz, nE ) + 1 )
          allocate ( work( lwork ) )
          ifail = 0
          call db2ink( Lzv1, nLz, Eval, nE, period, nLz, 4, 4, tL, tE,
     +                 Bc2, work, ifail )
          if( ifail .gt. 1 )then
            if( master )print *, 'iflag =', ifail
            call crash( 'TAUTAB', 'DB2INK failed' )
          end if
          if( master )print *, 'TAUTAB: table now ready'
        end if
        iuse = icmp
      end if
c Maclaurin/Freeman/Kalnajs disc
      if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'MFK ' ) )then
        tautab = sqrt( 4. * pi / 3. )
c uniform sphere
      else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'UNIS' ) )then
        tautab = pi
c Mestel/Toomre/Zang disc
      else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'MTZ ' ) )then
        Lz1 = Lz * Lzm * exp( .5 - E )
        tautab = algrng2( Eval, period1, nE, Lz1 )
        tautab = tautab * exp( E - 1. )
c isochrone disc
      else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'ISOC' ) )then
        tautab = algrng2( Eval, period1, nE, E )
c power law disc
      else if( ( ncmp .eq. 1 ) .and. ( ctype( 1 ) .eq. 'POWE' ) )then
        Lz1 = Lz * Lzm / Lzmax( E )
        tautab = algrng2( Eval, period1, nE, Lz1 )
        r = ( 1.d0 - dfcns( 3, icmp ) ) / ( 2.d0 * dfcns( 3, icmp ) )
        tautab = tautab * abs( E )**r
c general potential
      else
        if( E .ne. Elast )then
          Lzl = Lzmin( E )
          r = rcofe( E )
          Lzm = r * vcirc( r )
          Lzm = max( Lzm, Lzl + 1.d-8 )
          ak = akappa( r )
          Elast = E
        end if
c particle at rest in the centre
        if( Lzm .eq. 0. )then
          tautab = period( 1, 1 ) / ak
        else
c interpolate value from table
          Lz1 = ( abs( Lz ) - Lzl ) / ( Lzm - Lzl )
          E1 = max( E, Emine )
          E1 = min( E1, Emaxe )
          tautab = db2val(
     +                 Lz1, E1, 0, 0, tL, tE, nLz, nE, 4, 4, Bc2, work )
          tautab = tautab / ak
        end if
      end if
      return
      end
