      module dfarrays
c Copyright (C) 2015, Jerry Sellwood
      real*8, allocatable, save :: RShu( : ), svR( : )
      real*8, allocatable, save :: svRc( : ), svRlam( : )
      real*8, allocatable, save :: Ltab( : ), lhs( : ), rhs( : )
      real*8, allocatable, save :: Flam( : ), Fc( : )
      real*8, allocatable, save :: absc( :,: ), wgt( :,: ), Shukrn(:,: )
      real*8 norm
      integer npts, nrads
      parameter ( npts = 64, nrads = 200 )
      end module dfarrays

      real*8 function dfShu( E, Lz )
      use dfarrays
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c Distribution function of the form F_0(Lz) exp[-(E-E_c)/sigmu_u^2]
c   as Shu (1969, ApJ v158, p505)
c need to pre-tabulate the radial velocity dispersion in order to avoid
c   recursive calls to external routines
c   uses CMLIB routines DB2INK and DB2VAL
c   uses Burkardt's routine smsno
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
      real*8 upsilon
      common / dfsetup / upsilon
c
c externals
      external functn, parmuf
      real*8 akappa, db2val, errfn, GSigma, Lzmin, Omegac, Phitot, rcirc
      real*8 rcofe, splint2, taper, vcirc
c
c local arrays
      real*8, allocatable :: Ftab( : ), d( : )
      integer, allocatable :: iw( : )
      real*8, allocatable :: work( : )
c tables for 2D spline
      integer lwork, nE, nLz
      parameter ( nE = 101, nLz = 101 )
      real*8, allocatable :: Eval( : ), Lzv1( : ), table( :, : )
      real*8, allocatable :: Bc2( : ), tL( : ), tE( : )
      save Bc2, tE, tL, work
c
c local variables
      integer i, ia, iE, ir, ifail, iLz, iuse, iunit, ival, liw, lw
      integer nf, nval, uiparm
      logical lmake
      real*8 ak, arg, dR, E1, Ec, Elast, Elim, Ep, F, L, Ls, Lz1, Lzl
      real*8 Lzm, Om, Phi, Phieff, Q, R, res2, Rg, svRR, x, urparm, zero
      save Elast, Elim, iuse, Lzl, Lzm
      include 'inc/pi.f'
      parameter ( nval = 21, zero = 0 )
      data iuse / 0 /, Elast / -100. /
c
      if( iuse .ne. icmp )then
c create a 2D table of values if they are not aleady in a .dft file
        iunit = -1
        call opnfil( iunit, 'dft', 'unformatted', 'new', 'seq', i )
        lmake = i .eq. 0
        if( lmake )then
          if( master )print *,
     +            'No existing .dft file, starting set up of a Shu DF'
          allocate ( RShu( nrads ) )
          allocate ( svR( nrads ) )
c flag this as a live mass component
          rigidp( icmp ) = .false.
c ensure outer taper has at least a small radial range
          Lztmx( icmp ) = rtrunc( icmp )
          Lztmn( icmp ) = min( Lztmn( icmp ), .95 * Lztmx( icmp ) )
          if( master )print *, 'Disc tapered in radius from',
     +               sngl( Lztmn( icmp ) ), ' to', sngl( Lztmx( icmp ) )
c suppress any retrograde stars for now
          Ls = Lzcrt( icmp )
          Lzcrt( icmp ) = 0
c set table of svR values
          dR = rtrunc( icmp ) / real( nrads - 1 )
          Q = dfcns( 1, icmp )
          do ir = 1, nrads
            R = dR * real( ir - 1 )
            Om = Omegac( R )
            ak = akappa( R )
            ak = min( ak, 2. * Om )
            RShu( ir ) = R
            svR( ir ) = 3.36 * Q * GSigma( R ) * taper( R ) / ak
          end do
c extrapolate to R = 0 for exp disk only - funny values of kappa can creep in
          if( ctype( icmp ) .eq. 'EXP ' .and. ( .not. numcal ) )then
            x = svR( 2 ) - RShu( 2 ) * ( svR( 3 ) - svR( 2 ) ) /
     +                                 ( RShu( 3 ) - RShu( 2 ) )
            svR( 1 ) = max( svR( 1 ), x )
          end if
c initialize spline interpolation
          allocate ( svRc( nrads + 4 ) )
          allocate ( svRlam( nrads + 4 ) )
          R = RShu( nrads / 2 )
          svRR = splint2( RShu, svR, nrads, R, svRlam, svRc, .true. )
c set up rhs
          allocate ( rhs( nrads ) )
          do ir = 2, nrads
            R = RShu( ir )
            rhs( ir ) = GSigma( R ) * taper( R ) * R
          end do
c revise outer taper rule from radii to Lz
          r = rtrunc( icmp )
          L = r * vcirc( r )
          Lztmx( icmp ) = L
          r = Lztmn( icmp )
          Lztmn( icmp ) = r * vcirc( r )
c initial guess
          allocate ( Ftab( nval ) )
          allocate ( d( nval ) )
          allocate ( Ltab( nval ) )
          do ival = 2, nval
c linear spacing in L
            L = real( ival - 1 ) / real( nval - 1 )
            L = .5 * Lztmx( icmp ) * ( 1.d0 - cos( pi * L ) )
            Ltab( ival ) = L
            Rg = rcirc( L )
            Ftab( ival ) = 1
            d( ival ) = 1
          end do
c central values
          R = 0
          rhs( 1 ) = gsigma( R )
          svRR = splint2( RShu, svR, nrads, R, svRlam, svRc, .false. )
          x = ( Phitot( rmax ) - Phitot( 0.d0 ) ) / svRR**2
          norm = pi * svRR**2 * ( 1.d0 - exp( -x ) )
          Ltab( 1 ) = 0
          Ftab( 1 ) = rhs( 1 ) / norm
          d( 1 ) = 1
c entropy penalty parameter
          upsilon = 5.d-2        ! value that works best for the exp disk
c make table of Shukrn
          allocate ( absc( npts, nrads ) )
          allocate ( wgt( npts, nrads ) )
          allocate ( Shukrn( npts, nrads ) )
          do ir = 2, nrads
            R = RShu( ir )
c compute largest L whose orbit crosses R
            Phi = Phitot( R ) - Emaxe
            if( R .le. 0.d0 )then
              L = 0
            else if( R .ge. rmax )then
              L = Lztmx( icmp )
            else
              L = Phi / ( 1. / rmax**2 - 1. / R**2 )
              L = sqrt( 2. * L )
              L = min( L, Lztmx( icmp ) )
            end if
c set integration weights for Gauss-Legendre quadrature
            call GLqtab( zero, L, npts, wgt( 1, ir ), absc( 1, ir ) )
c integrate over L
            do ia = 1, npts
              L = absc( ia, ir )
              Rg = rcirc( L )
              arg = max( Rg, RShu( 1 ) )
              svRR = 0
              if( arg .le. RShu( nrads ) )svRR =
     +           splint2( RShu, svR, nrads, arg, svRlam, svRc, .false. )
              Shukrn( ia, ir ) = 0
              if( svRR .gt. 0.d0 )then
c argument of error fn = sqrt(-Phi_eff/sigma^2)
                Phieff = .5 * ( L / R )**2 + Phi
                arg = max( -Phieff, 0.d0 )
                arg = sqrt( arg ) / svRR
                ifail = 0
                Shukrn( ia, ir ) = errfn( arg, ifail )
c exponential factor
                Ec = .5 * ( L / Rg )**2 + Phitot( Rg ) - Emaxe
                Ep = Phieff - Ec
                arg = max( Ep, 0.d0 ) / svRR**2
                Shukrn( ia, ir ) = Shukrn( ia, ir ) * exp( -arg )
c remaining factors
                Shukrn( ia, ir ) =
     +                         sqrt( 2. * pi ) * Shukrn( ia, ir ) * svRR
              end if
            end do
          end do
c values for smsno
          liw = nval + 2
          lw = ( nval * ( nval - 1 ) ) / 2 + 12 * nval
          liw = max( liw, 60 )
          lw = max( lw, 77 + nval * ( nval + 17 ) / 2 )
          allocate ( iw( liw ) )
          allocate ( work( lw ) )
          allocate ( lhs( nrads ) )
c initialize
          call deflt( 2, iw, liw, lw, work )
c flag this as a cold start
          iw( 1 ) = 12
c allow more iterations
          iw( 18 ) = nval * 500
c parameters to suppress all output
          do i = 19, 24
            iw( i ) = 0
          end do
          iw( 21 ) = 0
          iw( 23 ) = -1
c least squares fitting routine
          call smsno( nval, d, Ftab, functn, iw, liw, lw, work,
     +                         uiparm, urparm, parmuf )
          if( master )print *, 'DFSHU: return status for smsno', iw( 1 )
          call functn( nval, Ftab, nf, res2, uiparm, urparm, parmuf )
c restore retrograde stars
          Lzcrt( icmp ) = Ls
c set disk energy bound and finish set up
          r = rtrunc( icmp )
          Elim = Phitot( R ) + .5 * vcirc( r )**2
          iuse = icmp
          if( master )print *,
     +                       'DFSHU: set up completed, building a table'
        end if
c allocate space
        allocate ( Eval( nE ) )
        allocate ( Lzv1( nLz ) )
        allocate ( table( nLz, nE ) )
        if( lmake )then
c choose Lzvals - concentrated towards the ends of the range
          do iLz = 1, nLz
            Lz1 = real( iLz - 1 ) / real( nLz - 1 )
            Lz1 = .5 * ( 1.d0 - cos( pi * Lz1 ) )
            Lz1 = max( Lz1, 0.d0 )
            Lz1 = min( Lz1, 1.d0 )
            Lzv1( iLz ) = Lz1
          end do
c choose Evals
          do iE = 1, nE
            E1 = real( iE - 1 ) / real( nE - 1 )
            E1 = .5 * ( Elim + Emine -
     +                ( Elim - Emine ) * cos( pi * E1 ) )
            E1 = max( E1, Emine )
            E1 = min( E1, Elim )
            Eval( iE ) = E1
            R = rcofe( E1 )
            Lzm = R * vcirc( R )
            Lzl = Lzmin( E1 )
            Lzl = min( Lzm, Lzl )
            Lzm = max( Lzm, Lzl + 1.d-8 )
c choose Lvals 
            do iLz = 1, nLz
              Lz1 = Lzv1( iLz )
              Lz1 = Lz1 * Lzm + ( 1. - Lz1 ) * Lzl
              Lz1 = max( Lz1, Lzl )
              Lz1 = min( Lz1, Lzm )
c find home radius and circular energy
              Rg = rcirc( Lz1 )
              Rg = min( Rg, RShu( nrads ) )
              Phi = Phitot( Rg )
              Ec = Phi
              if( Rg .gt. 0.d0 )Ec = Phi + 0.5 * ( Lz1 / Rg )**2
              Ec = min( Ec, E1 )
c look up required vely dispersion
              x = max( Rg, RShu( 1 ) )
              svRR = 0
              if( x .lt. RShu( nrads ) )svRR = splint2(
     +                      RShu, svR, nrads, x, svRlam, svRc, .false. )
              dfShu = 0
              if( svRR .gt. 0.d0 )then
c interpolate
                F = splint2( Ltab, Ftab, nval, Lz1, Flam, Fc, .false. )
c energy term
                dfShu = F * exp( ( Ec - E1 ) / svRR**2 )
c tabulate function for a full-mass disk for consistency with other DFs
                dfShu = dfShu / cmpmas( icmp )
              end if
c save this value
              table( iLz, iE ) = dfShu
            end do
          end do
          print *, 'DFSHU: table ready'
          write( iunit )table, Lzv1, Eval, Elim
        else
          call opnfil( iunit, 'dft', 'unformatted', 'old', 'seq', i )
          read( iunit )table, Lzv1, Eval, Elim
        end if
        close( iunit )
c fit a 2D cubic spline
        allocate ( tL( nLz + 4 ) )
        allocate ( tE( nE + 4 ) )
        allocate ( Bc2( nLz * nE ) )
        if( allocated( work ) )deallocate ( work )
        lwork = nLz * nE + 8 * ( max( nLz, nE ) + 1 )
        allocate ( work( lwork ) )
        ifail = 0
        call db2ink( Lzv1, nLz, Eval, nE, table, nLz, 4, 4, tL, tE,
     +               Bc2, work, ifail )
        if( ifail .gt. 1 )then
          if( master )print *, 'iflag =', ifail
          call crash( 'DFSHU', 'DB2INK failed' )
        end if
        if( master )print *, 'DFSHU: ready to use'
        iuse = icmp
      end if
c default value
      dfShu = 0
      if( E .lt. Elim )then
        if( E .ne. Elast )then
          Lzl = Lzmin( E )
          r = rcofe( E )
          Lzm = r * vcirc( r )
          Lzm = max( Lzm, Lzl + 1.d-8 )
          Elast = E
        end if
c star at rest in the center
        if( Lzm .eq. 0. )then
          dfShu = table( 1, 1 )
        else
c interpolate value from table
          E1 = max( E, Emine )
          E1 = min( E1, Elim )
          Lz1 = ( abs( Lz ) - Lzl ) / ( Lzm - Lzl )
          Lz1 = max( Lz1, 0.d0 )
          Lz1 = min( Lz1, 1.d0 )
          dfShu = db2val(
     +                 Lz1, E1, 0, 0, tL, tE, nLz, nE, 4, 4, Bc2, work )
        end if
        dfShu = max( dfShu, 0.d0 )
c apply retrograde star rule
        if( Lz .lt. -Lzcrt( icmp ) )then
          x = 0.
        else if( Lz .lt. Lzcrt( icmp ) )then
          x = Lz / Lzcrt( icmp )
          x = .5 + .75 * x - .25 * x**3
        else
          x = 1
        end if
        dfShu = x * dfShu
      end if
      return
      end

      subroutine functn( nval, Ftab, nf, T, uiparm, urparm, ufparm )
      use dfarrays
c  Copyright (C) 2015, Jerry Sellwood
c  For use with Burkardt's smsno, which does not require gradients
      implicit none
c
c calling arguments
      integer nf, nval, uiparm
      external ufparm
      real*8 T, urparm, Ftab( nval )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      real*8 upsilon
      common / dfsetup / upsilon
c
c externals
      real*8 GSigma, splint2
c
c local variables
      integer ia, icall, ir, ival
      logical firstc
      real*8 alpha, C, F, L, minT, R, S, zero
      parameter ( alpha = -1.d+5, zero = 0 )
      include 'inc/pi.f'
      save firstc, icall, minT
c
      data firstc, icall, minT / .true., 0, 1.d+5 /
c
      if( firstc )then
        allocate ( Flam( nval + 4 ) )
        allocate ( Fc( nval + 4 ) )
        firstc = .false.
      end if
c fit an interpolating spline to the new table
      L = Ltab( nval / 2 )
      F = splint2( Ltab, Ftab, nval, L, Flam, Fc, .true. )
c compute least squares difference from predicted surface density at each R
      C = 0
      do ir = 2, nrads
        lhs( ir ) = 0
        do ia = 1, npts
          L = absc( ia, ir )
          F = splint2( Ltab, Ftab, nval, L, Flam, Fc, .false. )
          lhs( ir ) = lhs( ir ) + wgt( ia, ir ) * F * Shukrn( ia, ir )
        end do
c sum squares of relative error - ignoring outer taper
        R = RShu( ir )
        C = C + ( ( lhs( ir ) - rhs( ir ) ) / ( R * GSigma( R ) ) )**2
      end do
      C = C + ( ( Ftab( 1 ) * norm - rhs( 1 ) ) / rhs( 1 ) )**2
c compute negative of entropy term
      S = 0
      do ival = 1, nval
        F = Ftab( ival )
        if( F .gt. 0.d0 )then
          S = S + F * log( F )
        else
c add a large penalty for a negative value
          S = S + alpha * F
        end if
      end do
c combine sum of squares and penalty term
      T = C + upsilon * S / real( nval )
      minT = min( T, minT )
      icall = icall + 1
      if( ( mod( icall, 500 ) .eq. 1 ) .and.
     +     master )print *, icall, 'best so far', minT
      return
      end

      subroutine parmuf
      print *, 'parmuf called!'
      return
      end
