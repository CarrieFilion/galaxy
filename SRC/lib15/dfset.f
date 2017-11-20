      subroutine dfset
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c main routine for deterministic generation of particles from a DF that can
c   be expressed as a function of energy and angular momentum
c The routine selects energy values between emin and emax evenly spaced
c   in mass.  For each of these energies, E, it then selects values of
c   angular momentum between Lzmin(E) & Lzmax(E) also equally spaced in
c   in mass, but with the first fraction selected randomly.
c Each of the selected (E,h) values defines an orbit in the potential and the
c   phases of particles on this orbit are selected by subroutine ORBITN
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
      include 'inc/orbval.f'
c
c externals
      logical gtlogl
      real*8 dblint, dfwght, gofx, Lzmax, Lzmin, psranu, ranuni, satab
      real*8 sglint, tautab
c
c local arrays
      integer jdim
      parameter ( jdim = 500 )
      real star( jdim )
      real*8 w( 4 )
c
c local variables
      integer i, iE, ifail, iLz, ind, irec, j, k, ner, nrec, nw
      logical check
      real rnew
      real*8 actmas, E, El, Eu, E1, E2, f, fract, Lz, Lzbase, LzEmax
      real*8 LzEmin, Lzf, Lzmas, Lz1, Lz2, tol
c
c initialize tables
      Lz = Lzmin( Emine )
      if( cdft( icmp ) .eq. 'KALN' )f = gofx( Lz )
      if( sphrod( icmp ) )then
        f = satab( Emine, Lz )
      else
        f = tautab( Emine, Lz )
      end if
      if( cmprssd )then
        if( .not. disc( icmp ) )then
c set icomp
          do i = 1, ncomp
            if( iccmp( i ) .eq. icmp )icomp = i
          end do
c enlarge outer radius as needed
          rnew = 1000
          f = 0
          if( r1cusp )f = -10
          j = iccmp( icomp )
          do i = nradad, 1, -1
            if( rhon( i, icomp ) .le. f )rnew = arad( i )
          end do
          if( ( ncmp .gt. 1 ) .and. ( rnew .gt. rtrunc( icmp ) ) )then
           print *, 'rtrunc increased from', rtrunc( icmp ), ' to', rnew
c           rtrunc( icmp ) = rnew
           call cutoff( rnew )
          end if
        else
          icomp = 0
          rnew = rtrunc( icmp )
          call cutoff( rnew )
        end if
      end if
c find total active mass
      actmas = dblint( Emine, Emaxe )
      check = ( .not. uqmass )
      if( cmprssd .and. ( icomp .gt. 0 ) )then
        f = ( gmn( nradad, icomp ) - actmas ) / gmn( nradad, icomp )
        check = check .and. ( abs( f ) .gt. 1.e-3 )
        if( check )then
          print *, 'expected value', sngl( gmn( nradad, icomp ) )
          print *, 'integrated value', sngl( actmas )
          print *, 'percentage mass error =', 100. * sngl( f )
        end if
      else
c set fmass if not already set (eg from .pft file)
        if( fmass( icmp ) .ne. 0. )then
          tol = fmass( icmp ) - actmas / cmpmas( icmp )
          check = check .and. ( abs( tol ) .gt. 1.d-2 )
          if( check )print *, 'Change of active mass fraction',
     +       fmass( icmp ), sngl( actmas / cmpmas( icmp ) ), sngl( tol )
        end if
      end if
      if( check )then
        if( .not. gtlogl( 'Is this acceptable' )
     +            )call crash( 'DFSET', 'Intolerable mass discrepancy' )
      end if
      fmass( icmp ) = actmas / cmpmas( icmp )
      print
     +  '( ''0Active mass fraction of component'', i2, '' ='', f10.4 )',
     +                                               icmp, fmass( icmp )
c write header record
      if( sphrod( icmp ) )then
        nw = 4
      else
        nw = 3
      end if
      if( uqmass )nw = nw + 1
      irec = nw * ( jdim / nw )
      call dfhead( ndistf, irec, .true., .true. )
c
      k = 0
      nrec = 0
      ner = jdfcs( 1, icmp ) / 20
      ner = max( ner, 1 )
c choose values of E
      Eu = Emine
      do iE = 1, jdfcs( 1, icmp )
        El = Eu
        f = actmas * ( dble( iE ) ) / dble( jdfcs( 1, icmp ) )
        f = min( f, actmas )
c set limits
        E1 = Emine
        E2 = Emaxe
        fract = f
        w( 1 ) = f - actmas
        ind = -1
        tol = 1.e-6 * ( Emaxe - Emine )
        ifail = 1
c find zero of function
        do while ( ind .ne. 0 )
          call fndzro( E1, E2, fract, tol, 0, w, ind, ifail )
          if( ind .ne. 0 )fract = f - dblint( Emine, E1 )
        end do
        if( ( ifail .gt. 0 ) .and. ( ifail .lt. 4 ) )then
          call crash( 'DFSET', 'Failed to find zero of M(E)' )
        end if
        Eu = E1
        if( mod( ie, ner ) .eq. 0 )print
     +                '( '' Energy cut no'', i7, f10.5 )', ie, Eu
        Lzbase = ranuni( 0. )
c start at a random point in the cycle of pseudo random numbers
        j = Lzbase * real( jdfcs( 2, icmp ) ) + 1
        do i = 1, j
          f = psranu( jdfcs( 2, icmp ) )
        end do
c choose El < E < Eu at a random fraction of mass for each Lz
        do iLz = 1, jdfcs( 2, icmp )
          f = psranu( jdfcs( 2, icmp ) )
          f = dblint( Emine, El ) * ( 1. - f ) +
     +        dblint( Emine, Eu ) * f
c find zero of function
          E1 = El
          E2 = Eu
          ind = 1
          ifail = 1
          do while ( ind .ne. 0 )
            call fndzro( E1, E2, fract, tol, 0, w, ind, ifail )
            if( ind .ne. 0 )fract = f - dblint( Emine, E1 )
          end do
          if( ( ifail .gt. 0 ) .and. ( ifail .lt. 4 ) )then
            print *, 'ifail =', ifail
            print '( 3f12.6, 2e15.5 )', E1, E2, f,
     +           f - dblint( Emine, E1 ), f - dblint( Emine, E2 )
            call crash( 'DFSET', 'Failed to find zero of M(Lz) 1' )
          end if
          E = E1
c integrate distfn along cut at this E
          LzEmax = Lzmax( E )
          LzEmin = max( -LzEmax, Lzrmn( icmp ) )
          if( LzEmin .eq. 0. )LzEmin = Lzmin( E )
          Lzmas = sglint( LzEmin, LzEmax, E )
c choose value of Lz
          Lzf = Lzmas *
     +               ( dble( iLz ) - Lzbase ) / dble( jdfcs( 2, icmp ) )
c set limits
          Lz1 = LzEmin
          Lz2 = LzEmax
          fract = Lzf
          w( 1 ) = Lzf - Lzmas
          ind = -1
          tol = 1.e-6 * ( LzEmax - LzEmin )
          ifail = 1
c find zero of function
          do while ( ind .ne. 0 )
            call fndzro( Lz1, Lz2, fract, tol, 0, w, ind, ifail )
            if( ind .ne. 0 )fract = Lzf - sglint( LzEmin, Lz1, E )
          end do
          if( ( ifail .gt. 0 ) .and. ( ifail .lt. 4 ) )then
            print *, Lz1, Lz2, fract
            call crash( 'DFSET', 'Failed to find zero of M(Lz) 2' )
          end if
          Lz = Lz1
c choose sets of coordinates for this orbit
          do j = 1, jdfcs( 3, icmp )
c            if( ie .lt. 0 )then
c              ri = ranuni( 0. )
c              ui = ranuni( 0. )
c            else
              call orbitn( E, Lz )
c            end if
            star( k + 1 ) = ri
            star( k + 2 ) = ui
            star( k + 3 ) = vi
            if( sphrod( icmp ) )then
              star( k + 4 ) = zi
              k = k + 4
            else
              k = k + 3
            end if
            if( uqmass )then
              star( k + 1 ) = fmass( icmp ) / dfwght( E, Lz )
              k = k + 1
            end if
            if( k + nw .gt. irec )then
              write( ndistf )( star( i ), i = 1, k )
              k = 0
            end if
          end do
        end do
      end do
c write out remaining data
      if( k .gt. 0 )write( ndistf )( star( i ), i = 1, k )
      return
      end
