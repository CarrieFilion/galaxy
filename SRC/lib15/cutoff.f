      subroutine cutoff( rcut )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c sets up common blocks elims and taper to restrict the energy and angular
c   momentum of particles to remain within the specified cut off radius
c
c calling argument
      real rcut
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
      real*8 Phitot, vcirc
c
c local variables
      integer ip, jcmp, jp, kp
      logical lknwn
      real*8 E, r
      include 'inc/pi.f'
c
      jcmp = icmp
      lknwn = .true.
      if( snglph )then
        jp = icmp
        kp = icmp
      else
        jp = 1
        kp = ncmp
      end if
c work over components
      do ip = jp, kp
        if( disc( ip ) )then
c trap finite discs
          if( imtyp( ip ) .gt. 22 )then
            call crash( 'CUTOFF', 'Unknown DTYPE' )
          else if( ( ctype( ip ) .eq. 'MFK ' ) .or.
     +             ( ctype( ip ) .eq. 'HUNT' ) .or.
     +             ( ctype( ip ) .eq. 'COSI' ) .or.
     +             ( ctype( ip ) .eq. 'FMES' ) .or.
     +             ( ctype( ip ) .eq. 'SPLY' ) .or.
     +             ( ctype( ip ) .eq. 'PPLY' ) .or.
     +             ( ctype( ip ) .eq. 'ARBT' ) .or.
     +             ( ctype( ip ) .eq. 'UNKN' ) )then
            trnctd( ip ) = .false.
          else
            trnctd( ip ) = rcut .lt. 999.
            if( master )then
              if( trnctd( ip ) )then
                if( rcut .ne. sngl( rmax ) )write( no,
     +                   '( '' Truncation radius set to'', f8.2 )' )rcut
              else
                if( sngl( rmax ) .lt. 999. )write( no,
     +                         '( '' Outer truncation eliminated'' )' )
              end if
            end if
          end if
        end if
c determine whether energy bounds can be set
        lknwn = lknwn .and. ( ctype( ip ) .ne. 'UNKN' )
      end do
c set up energy limits etc if possible
      if( lknwn )then
        rmax = rcut
        Phimax = Phitot( rmax )
        Emaxe = Phimax
c maximum energy for a circular orbit at the edge of a disk
        do ip = jp, kp
          if( disc( ip ) )then
            r = min( rtrunc( ip ), rcut )
            E = Phitot( r ) + .5 * vcirc( r )**2
            Emaxe = max( Emaxe, E )
          end if
        end do
c energy of a circular orbit at the edge of the central hole
        Emine = Phitot( rhole ) + .5 * vcirc( rhole )**2
c trap nonsense
        if( Emine .gt. Emaxe )then
          print *, jcmp, ncmp, sngl( Emine ), sngl( Emaxe )
          call crash( 'CUTOFF', 'Emine > Emaxe' )
        end if
c work over mass components
        do icmp = jp, kp
          if( dist( icmp ) )then
c Emax is less for Omega models
            if( cdft( icmp ) .eq. 'OMEG' )Emaxe = Phi0 +
     +         .5 * ( ( 3. * pi / 4. ) + sqrt( 3. * pi / 4. ) * omegak )
c Emax is set in .pft header for DFIT models
            if( ( cdft( icmp ) .eq. 'DFIT' ) .or.
     +          ( cdft( icmp ) .eq. 'ADIA' ) )Emaxe = dfcns( 3, icmp )
c Emax for the Polyachenko-Shukhman DF for a uniform sphere
            if( cdft( icmp ) .eq. 'USPS' )Emaxe =
     +                                     Emaxe + .5 * vcirc( rmax )**2
c set up default tapers (if not already set)
            if( .not. Lztapr( icmp ) )then
              Lztmx( icmp ) = rmax * vcirc( rmax )
              Lztmn( icmp ) = Lztmx( icmp )
              if( ( Lzrmn( icmp ) .ge. 0.d0 ) .and. dist( icmp ) )then
c most retrograde value
c
c Kalnajs fn with Agris's retrograde rule
                if( ( cdft( icmp ) .eq. 'KALN' ) .and. retract )then
                  Lzrmn( icmp ) = -Lztmx( icmp )
c Kalnajs, Lia, Evans & Collett or Sawamura functions
                else if( ( cdft( icmp ) .eq. 'KALN' ) .or.
     +                   ( cdft( icmp ) .eq. 'LIA ' ) .or.
     +                   ( cdft( icmp ) .eq. 'EVCO' ) .or.
     +                   ( cdft( icmp ) .eq. 'SAWA' ) )then
                  Lzrmn( icmp ) = -min( Lzcrt( icmp ), Lztmx( icmp ) )
c Miyamoto function
                else if( cdft( icmp ) .eq. 'MIYA' )then
                  Lzrmn( icmp ) = -Lztmx( icmp )
c Omega model
                else if( cdft( icmp ) .eq. 'OMEG' )then
                  Lzrmn( icmp ) =
     +             .5 * ( dfcns( 3, icmp ) - 1. ) * sqrt( 3. * pi / 4. )
c no retrograde stars for power law discs or spherical models
                else if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +                   ( cdft( icmp ) .eq. 'EVAN' ) .or.
     +                   ( cdft( icmp ) .eq. 'DJDZ' ) .or.
     +                   ( cdft( icmp ) .eq. 'DFIT' ) .or.
     +                   ( .not. sphrod( icmp ) ) )then
                  Lzrmn( icmp ) = 0
                else
                  call crash( 'CUTOFF', 'Unrecognized dftype' )
                end if
              end if
            end if
          end if
        end do
      end if
c restore icmp
      icmp = jcmp
      return
      end
