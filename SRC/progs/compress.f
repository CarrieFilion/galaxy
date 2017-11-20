      program compress
      use aarrays
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c program to implement Young's (1980) algorithm for compression of a
c   spherical system with a known DF by the addition of rigid matter
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
      external lrhohal, lPhi, gmassh, vcirc
      logical gtlogl
      real*8 distfn, E0tab, gmassh, gmext, Phihal, rhohal
c
c local variables
      integer i, j
      logical convrg, conserve, cuspy
      real errm( 2 ), x
      real*8 diffl, diffm, E, errmas, r, rm, tm, totmas( 2 ), extmas
      parameter ( conserve = .false. )
c
c set defaults
      call setnag
      no = 6
      master = .true.
      call boilrp( .true. )
      call gtintg( 'Enter 0 for terminal input', ni )
      if( ni .ne. 0 )call getset
c set initial model
      call msetup
      nradad = 1001
c determine type of compress
      ncomp = 0
      nrigid = 0
      extmas = 0
      cuspy = .false.
      r2cusp = .false.
      snglph = .true.
      i = 0
      do icmp = 1, ncmp
        if( idftyp( icmp ) .gt. 1 )then
          rigidp( icmp ) = .false.
          ncomp = ncomp + 1
          if( ncomp .gt. 2
     +       )call crash( 'COMPRESS', 'More than 2 compressible halos' )
c remember properties of uncompressed halo(s)
          iccmp( ncomp ) = icmp
          j = imtyp( icmp )
          ihalo0( ncomp ) = j
          cuspy = cuspy .and. ( ( j .eq. 4 ) .or. ( j .eq. 6 ) .or.
     +                        ( j .eq. 16 ) )
          cuspy = cuspy .or. ( ( j .eq. 9 ) .or. ( j .eq. 20 ) .or.
     +                       ( j .eq. 25 ) ) .or. ( j .eq. 32 )
          if( j .eq. 9 )r2cusp = .true.
          idfn0( ncomp ) = idftyp( icmp )
          dfm0( ncomp ) = dfcns( 1, icmp )
          dfb0( ncomp ) = dfcns( 2, icmp )
          cc0( ncomp ) = dfcns( 3, icmp )
          r = rtrunc( icmp )
          Emaxo( ncomp ) = Phihal( r )
        else
          nrigid = nrigid + 1
          rigidp( icmp ) = .true.
          if( nrigid .gt. 2
     +         )call crash( 'COMPRESS', 'More than 2 rigid components' )
c remember compressing mass component(s)
          ircmp( nrigid ) = icmp
          extmas = extmas + cmpmas( icmp )
c set disk thickness and functional form - Gaussian with z0=rscale/10
          if( disc( icmp ) )then
            i = i + 1
            if( i .gt. 1 )call crash( 'COMPRESS', 'more than 1 disk' )
            z0init( icmp ) = .1 * rscale( icmp )
            iztyp( icmp ) = 2
            print *, 'Warning: disk thickness 0.1Rd and Gaussian' //
     +                                       ' profile are hard-wired'
            if( .not. gtlogl(
     +          'is this OK' ) )call crash( 'COMPRESS', 'user abort' )
c set up table for gmext
            r = .5 * rtrunc( icmp )
            x = gmext( r )
          end if
        end if
      end do
      if( ncomp .eq. 0 )call crash( 'COMPRESS', 'No compressible halo' )
      if( extmas .eq. 0.d0 )then
        print *, 'no rigid mass'
        if( .not. gtlogl( 'Do you want to continue' )
     +                        )call crash( 'COMPRESS', 'No rigid mass' )
      end if
c find max rad and select basis for radial spacing
      rm = 0
      do icomp = 1, ncomp
        icmp = iccmp( icomp )
        r = rtrunc( icmp )
        if( ihalo0( icomp ) .eq. 6 )r = 1.5 * r
        rm = max( r, rm )
      end do
      r1cusp = cuspy
      r0core = .not. cuspy
c allocate arrays
      allocate ( arad( mradad ) )
      allocate ( phin( mradad ) )
      allocate ( gmn( mradad, ncomp ) )
      allocate ( rhon( mradad, ncomp ) )
      allocate ( plamda( mradad + 4 ) )
      allocate ( phinc( mradad + 4 ) )
      allocate ( glamda( mradad + 4, ncomp ) )
      allocate ( gmnc( mradad + 4, ncomp ) )
      allocate ( rlamda( mradad + 4, ncomp ) )
      allocate ( rhonc( mradad + 4, ncomp ) )
c initialize tables
      do icomp = 1, ncomp
        icmp = iccmp( icomp )
c set revised energy bounds for this component only
        x = rtrunc( icmp )
        call cutoff( x )
c build tables of E_0(J1,J2) and df(E, Lz) for this component
        E = E0tab( 0.d0, 0.d0 )
        E = max( E, -1.d0 )
        E = distfn( E, 0.d0 )
c initialize halo density
        do i = 1, nradad
          if( icomp .eq. 1 )then
            r = dble( i - 1 ) / dble( 2 * ( nradad - 1 ) )
            if( r0core )then
              r = rm * r / ( 1. - r )
              r = max( r, 5.d-4 )
            else if( r1cusp .or. r2cusp )then
              r = rm * ( r / ( 1. - r ) )**2
              r = max( r, 1.d-7 )
            else
              call crash( 'COMPRESS', 'Unknown rule for abscissae' )
            end if
            arad( i ) = r
            Phin( i ) = Phihal( r )
          else
            r = arad( i )
          end if
          gmn( i, icomp ) = gmassh( r )
          rhon( i, icomp ) = rhohal( r )
          if( r1cusp .or. r2cusp
     +                   )rhon( i, icomp ) = log10( rhon( i, icomp ) )
        end do
        totmas( icomp ) = gmn( nradad, icomp )
        print *,
     +      'totmas for component', icmp, ' is', sngl( totmas( icomp ) )
      end do
c
c compute truncated mass
      call renew_rho( diffm )
c change halo type
      do i = 1, ncomp
        icmp = iccmp( i )
        imtyp( icmp ) = 23
        ctype( icmp ) = 'ADIA'
        idftyp( icmp ) = 25
        cdft( icmp ) = 'COMP'
      end do
      snglph = .false.
      cmprssd = .true.
c compute initial potential
      diffl = Phimax
      call renew_Phi( diffm )
c reset energy bounds for multi-component model
      x = arad( nradad )
      call cutoff( x )
      do icomp = 1, ncomp
        print *, 'Notional and truncated total masses',
     +           sngl( totmas( icomp ) ), sngl( gmn( nradad, icomp ) )
        totmas( icomp ) = gmn( nradad, icomp )
        icmp = iccmp( icomp )
      end do
c this assumes that the halo is the first compressible component
      if( conserve )then
        call crash( 'COMPRESS', 'Code to conserve mass not ready' )
        cmpmas( ncomp ) = ( totmas( 1 ) - extmas ) / totmas( ncomp )
      end if
c
c iterate
      iterad = 1
      convrg = .false.
      do while ( .not. convrg )
        print *, 'starting iteration', iterad
        call renew_rho( diffm )
        call renew_Phi( diffm )
c
        do i = 1, ncomp
          if( ( i .eq. 1 ) .and. conserve )then
            errmas = gmn( nradad, i ) - extmas
          else
            errmas = gmn( nradad, i )
          end if
          errm( i ) = ( totmas( i ) - errmas ) / totmas( i )
          print *,
     +      'Fractional error in mass of component', i, ' is', errm( i )
        end do
        convrg = diffm .lt. 2.d-5
        print *, 'Iteration', iterad, ' Max pot diff', sngl( diffm )
        iterad = iterad + 1
        if( ( iterad .gt. 2 ) .and. ( .not. convrg ) )then
          if( ( diffl - diffm ) / diffm .lt. .01 )then
            print *, 'The iteration is not making good progress'
            if( gtlogl( 'Do you want it to continue' ) )then
              continue
            else
              if( gtlogl( 'Is the current difference acceptable' ) )then
                convrg = .true.
              else
                call crash( 'COMPRESS', 'Run abandoned' )
              end if
            end if
          end if
        end if
        diffl = diffm
      end do
c check total mass
      tm = 0
      x = 0
      do i = 1, ncomp
        if( ( i .eq. 1 ) .and. conserve )then
          errmas = gmn( nradad, i ) - extmas
        else
          errmas = gmn( nradad, i )
        end if
        tm = tm + totmas( i )
        x = x + errmas
      end do
      if( conserve )errmas = errmas + extmas
      print *, 'Original and current total masses', sngl( tm ), x
      call compin( .false. )
      end

      real*8 function lPhi( a )
      implicit none
      real*8 a, r, Phihal
      r = 10.**a
      lPhi = -Phihal( r )
      if( lPhi .le. 0.d0 )then
        lPhi = -10
      else
        lPhi = log10( lPhi )
      end if
      return
      end

      real*8 function ldfnow( E )
      real*8 E, distfn, L, Lzmax
      L = .5 * Lzmax( E )
      ldfnow = distfn( E, L )
      if( ldfnow .gt. 0.d0 )then
        ldfnow = log10( ldfnow )
      else
        ldfnow = -10
      end if
      return
      end

      subroutine renew_Phi( difm )
      use aarrays
      implicit none
c
c calling argument
      real*8 difm
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gmnx, Phinx, splint2
c
c local variables
      integer i
      real*8 r, s
c
      do icomp = 1, ncomp
        icmp = iccmp( icomp )
c compute new mass array
        do i = 1, nradad
          r = arad( i )
          gmn( i, icomp ) = gmnx( r )
        end do
c fit new spline
        r = arad( 1 )
        s = splint2( arad, gmn( 1, icomp ), nradad, r,
     +               glamda( 1, icomp ), gmnc( 1, icomp ), .true. )
      end do
c compute new potential array
      difm = 0
      rmax = arad( nradad )
      do i = 1, nradad
        r = arad( i )
        s = Phinx( r )
        difm = max( difm, abs( s - Phin( i ) ) )
        Phin( i ) = s
      end do
c fit new spline
      r = arad( 1 )
      s = splint2( arad, Phin, nradad, r, plamda, Phinc, .true. )
      return
      end

      real*8 function gmnx( r )
      implicit none
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 rhohal
c
c local arrays
      integer npts
      parameter ( npts = 64 )
      real*8 absc( npts ), wght( npts )
c
c local variables
      integer i, j, n, nval( 15 )
      real*8 rp
      include 'inc/pi.f'
c
      data nval / 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24, 32, 48, 64 /
c
      gmnx = 0
      if( r .gt. 0.d0 )then
c Gauss-Legendre quadrature over interior mass
        n = npts
c fewer abscissae needed at small radii, but must be an allowed value
        if( r .lt. 1.d0 )then
          n = r * real( npts )
          n = max( n, 2 )
          n = min( n, npts )
          j = 2
          i = 1
          do while ( j .lt. n )
            i = i + 1
            j = nval( i )
          end do
          n = j
        end if
        call GLqtab( 0.d0, r, n, wght, absc )
        do i = 1, n
          rp = absc( i )
          gmnx = gmnx + wght( i ) * rp**2 * rhohal( rp )
        end do
        gmnx = 4. * pi * gmnx
      end if
      return
      end

      real*8 function Phinx( r )
      implicit none
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      external frtot
      real*8 gmassh, Phiext, quad_gnr
c
c local variables
      integer ier
      real*8 epsr, rp
c
      Phinx = 0
      if( r .lt. rmax )then
c integrate central attraction
        epsr = 1.d-30
        Phinx = quad_gnr( frtot, r, rmax, epsr, epsr, ier )
        if( ier .gt. 2 )then
          print *, 'ier =', ier, ' from quad_gnr'
          if( ier .gt. 5 )call crash( 'PHINX', 'QUADPACK error' )
        end if
      end if
c add the contributions to infinity assuming no additional mass
      rp = max( r, rmax )
      do icomp = 1, ncomp
        icmp = iccmp( icomp )
        Phinx = Phinx - gmassh( rmax ) / rp
      end do
c add potential of compressing mass
      do irigid = 1, nrigid
        icmp = ircmp( irigid )
        Phinx = Phinx + Phiext( rp )
      end do
      return
      end

      subroutine renew_rho( difm )
      use aarrays
      implicit none
c
c calling argument
      real*8 difm
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 rhoint, splint2
c
c local variables
      integer i
      real r1
      real*8 r, rho, rhol
c
c re-initialize ajrtab
      iusead = 0
c work over compressible components
      do icomp = 1, ncomp
        icmp = iccmp( icomp )
c determine whether this is the most extensive component
        r1 = 0
        do i = 1, ncmp
          r1 = max( r1, rtrunc( i ) )
        end do
        if( cmprssd .and. ( rtrunc( icmp ) .eq. r1 ) )then
c outer boundary may need to be larger
          r1 = arad( nradad )
        else
c current truncation radius is large enough
          r1 = rtrunc( icmp )
        end if
c reset energy bounds
        call cutoff( r1 )
c update local density
        difm = 0
        do i = 1, nradad - 1
          r = arad( i )
          if( r .lt. rmax )then
            rho = rhoint( r )
c compute max change
            rhol = rhon( i, icomp )
            if( r1cusp .or. r2cusp )then
              if( rhol .gt. -8.d0 )then
                rhol = 10.**rhol
              else
                rhol = 0
              end if
            end if
            difm = max( difm, abs( rhol - rho ) )

            if( rho .lt. 0.d0 )then
              print *, 'Negative halo density!'
              print *, icomp, i,
     +               sngl( arad( i ) ), sngl( rho ), sngl( rhol )
              if( rho .lt. -1.d-4 )read *, rhol
              rho = 0
            end if

          else
            rho = 0
          end if
c save new value
          rhon( i, icomp ) = rho
          if( r1cusp .or. r2cusp )then
            if( rho .gt. 0.d0 )then
              rhon( i, icomp ) = log10( rho )
            else
              rhon( i, icomp ) = -10
            end if
          end if
        end do
        if( r1cusp .or. r2cusp )then
          rhon( nradad, icomp ) = -10
        else
          rhon( nradad, icomp ) = 0
        end if
c fit new spline
        r = arad( 1 )
        rho = splint2( arad, rhon( 1, icomp ), nradad, r,
     +                 rlamda( 1, icomp ), rhonc( 1, icomp ), .true. )
      end do
      return
      end

      real*8 function rhoint( r )
      implicit none
c Performs a double integration over all allowed velocities
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      real*8 E, rr, Pot
      common / intrho / E, rr, Pot
c
c externals
      external rhoifn
      real*8 distfn, Phitot, quad_Pat
c
c local arrays
      integer npts
      parameter ( npts = 64 )
      real*8 absc( npts ), wght( npts )
c
c local variables
      integer i, ier, nlim
      real*8 epsr, f, Lmax, u
      include 'inc/pi.f'
c
c set integration range
      rr = r
      Pot = Phitot( rr )
      if( Pot .lt. Phimax )then
c Gauss-Legendre quadrature over energies
        call GLqtab( Pot, Phimax, npts, wght, absc )
        rhoint = 0
        do i = 1, npts
          E = absc( i )
          if( rr .gt. 0.d0 )then
c adaptive integration over L
            Lmax = rr * sqrt( 2. * ( E - Pot ) )
            nlim = 0
            epsr = 1.e-10
            f = quad_Pat( rhoifn, 0.d0, Lmax, epsr, ier )
            rhoint = rhoint + wght( i ) * f
          else
            u = 2 * ( E - Pot )
            if( u .lt. 0. )call crash( 'RHOINT', 'E out of range' )
            u = sqrt( u )
            rhoint = rhoint + wght( i ) * u * distfn( E, 0.d0 )
          end if
        end do
      else
        rhoint = 0
      end if
      rhoint = 4. * pi * rhoint
      return
      end

      real*8 function rhoifn( L )
      implicit none
c
c calling argument
      real*8 L
c
c common blocks
c
      real*8 E, rr, Pot
      common / intrho / E, rr, Pot
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 distfn
c
c local variable
      real*8 u
c
      u = 2. * ( E - Pot ) - ( L / rr )**2
      if( u .le. 0. )then
        if( u .lt. -1.d-5 )then
          print *, E, Pot, L, rr, u
          call crash( 'RHOIFN', 'Argument out of range' )
        end if
        rhoifn = 0
      else
        u = sqrt( u )
        rhoifn = L * distfn( E, L ) / ( rr**2 * u )
      end if
      return
      end

      real*8 function Phirel( lr )
      implicit none
c
      real*8 lr, r, Phitot
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      r = 10.**lr
      Phirel = Phitot( r ) - Emine
      return
      end
