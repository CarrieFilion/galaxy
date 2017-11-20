      real*8 function dblint( E1, E2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the mass in a range of energies from direct integration of a DF
c   which is assumed to be of the form F(E,L)
c In order to avoid expensive repeated integrations when searching for the
c   energy of a particular mass fraction, we tabulate the results from each
c   small piece of the integral and interpolate.
c Integration by Simpson's rule appears to be accurate enough, since the
c   integrand is well behaved (ie, smooth and non-oscillatory).
c
c calling arguments
      real*8 E1, E2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 Lzmax, Lzmin, phitot, sglint, splint2
c
c local arrays
      integer idim, iuse, nval
      parameter ( idim = 1001 )
      real*8 c( idim + 4 ), Emass( idim ), Etab( idim ), lam( idim + 4 )
      save c, Emass, Etab, iuse, lam, nval
c
c local variables
      integer iE, iE1, nm
      logical even
      real*8 E, Elim, fE, f1, f2, f3, Lz1, Lz2
      include 'inc/pi.f'
c
      data iuse / 0 /
c
      if( iuse .ne. icmp )then
c create table of mass as a function of E - lower resolution OK for spheroids
        nval = idim
c        if( .not. disc( icmp ) )nval = 51
        nm = 2 * ( nval - 1 )
        Elim = Emaxe
        if( ctype( icmp ) .eq. 'DFIT' )Elim = dfcns( 3, icmp )
        fE = ( Elim - Emine ) / real( nm )
c integration using Simpson's rule
        f1 = 0.
        Etab( 1 ) = Emine
        Emass( 1 ) = 0.
        even = .true.
        do iE = 1, nm
          even = .not. even
          E = Emine + fE * real( iE )
          E = min( E, Elim )
c range of Lz
          Lz2 = Lzmax( E )
          Lz1 = max( -Lz2, Lzrmn( icmp ) )
          if( Lz1 .eq. 0. )Lz1 = Lzmin( E )
c special limit for a composite Omega model
          if( ( cdft( icmp ) .eq. 'CPB0' ) .and.
     +        ( E - phitot( 0.d0 ) .gt. .375 * pi ) )Lz1 = 0
          if( even )then
            f3 = sglint( Lz1, Lz2, E )
            iE1 = iE / 2
            Etab( iE1 + 1 ) = E
c Simpson's rule
            Emass( iE1 + 1 ) = Emass( iE1 ) +
     +                                   fE * ( f1 + 4. * f2 + f3 ) / 3.
            f1 = f3
          else
            f2 = sglint( Lz1, Lz2, E )
          end if
        end do
c initialize spline
        E = Etab( nval / 2 )
        f1 = splint2( Etab, Emass, nval, E, lam, c, .true. )
c set flag and constants
        iuse = icmp
      end if
c look up mass fractions in table
      if( E2 .gt. Etab( nval ) )then
        f2 = Emass( nval )
      else
        f2 = splint2( Etab, Emass, nval, E2, lam, c, .false. )
      end if
      if( E1 .lt. Etab( 1 ) )then
        f1 = Emass( 1 )
      else
        f1 = splint2( Etab, Emass, nval, E1, lam, c, .false. )
      end if
      dblint = f2 - f1
      return
      end
