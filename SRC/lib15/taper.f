      real*8 function taper( Lz )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c Returns a value between 0 and 1 with which to multiply the DF or surface
c   density.  The rules applied are for outer and (possibly) inner tapers
c   and a retrograde star rule
c
c The outer taper rule is either a cubic function from Lztmn to Lztmx or a
c   Zang-style power law about Lzmno
c
c The only inner taper rule is a Zang-style power law about Lz = 1
c
c There are a number of retrograde star rules
c
c calling argument
      real*8 Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      integer k
      logical done
      real*8 a, ti, to, x
c
      if( ( icmp .lt. 1 ) .or. ( icmp .gt. ncmp ) )call crash( 'TAPER',
     +                                        'Nonsense value of icmp' )
      done = .false.
c power law tapers for Zang models
      if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +    ( cdft( icmp ) .eq. 'EVAN' ) )then
        if( indexi( icmp ) .gt. 0. )then
          ti = Lz**indexi( icmp )
          ti = ti / ( 1. + ti )
        else
          ti = 1
        end if
        to = Lz / Lzmno( icmp )
        to = 1. / ( 1. + to**indexo( icmp ) )
        taper = ti * to
        done = .true.
c Evans-Collett function for the Rybicki disc
      else if( cdft( icmp ) .eq. 'EVCO' )then
        to = Lz / Lzmno( icmp )
        to = 1. / ( 1. + to**indexo( icmp ) )
        taper = to
      else
c general outer edge taper if required - Lz is interpreted as R for no DF
        taper = 1.
        if( trnctd( icmp ) )then
c set function to zero beyond maximum Lz
          if( abs( Lz ) .gt. Lztmx( icmp ) )then
            taper = 0.
c cubic taper near outer edge
          else if( abs( Lz ) .gt. Lztmn( icmp ) )then
            x = 1. - ( abs( Lz ) - Lztmn( icmp ) ) /
     +                                 ( Lztmx( icmp ) - Lztmn( icmp ) )
            taper = -2. * x**3 + 3. * x**2
          end if
          done = .true.
        end if
      end if
c
c retrograde star rules
c
c no DF - no sensible rule possible
      done = done .and. ( .not. dist( icmp ) )
      done = done .or. retract
      if( .not. done )then
c non-rotating Kalnajs function
        if( ( cdft( icmp ) .eq. 'KALN' ) .and. ( .not. retract ) .and.
     + ( dfcns( 1, icmp ) .eq. 3 ) .and. ( Lzcrt( icmp ) .gt. 0. ) )then
          taper = .5 * taper
          done = .true.
c Evans-Collett functions for the Rybicki disc
        else if( cdft( icmp ) .eq. 'EVCO' )then
c non-rotating case
          if( dfcns( 1, icmp ) .lt. 0. )then
            taper = .5 * taper
            done = .true.
c special case
          else if( ( Lz .eq. 0. ) .and. ( Lzcrt( icmp ) .eq. 0. ) )then
            taper = .5 * taper
            done = .true.
          end if
c no retrograde stars in power-law disks
        else if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +           ( cdft( icmp ) .eq. 'SHUE' ) .or.
     +           ( cdft( icmp ) .eq. 'EVAN' ) )then
          done = .true.
        end if
      end if
      if( .not. done )then
c cubic taper retrograde rule - Kalnajs, Lia, Dejonghe and Dejonghe-deZeeuw
        if( ( cdft( icmp ) .eq. 'KALN' ) .or.
     +      ( cdft( icmp ) .eq. 'LIA ' ) .or.
     +      ( cdft( icmp ) .eq. 'DEJO' ) .or.
     +      ( cdft( icmp ) .eq. 'DJDZ' ) )then
          if( Lz .lt. Lzcrt( icmp ) )then
            x = abs( Lz ) / Lzcrt( icmp )
            x = .5 - .75 * x + .25 * x**3
c x = .5 - 1.5 * x**2 + x**3
            if( Lz .le. 0. )taper = x * taper
            if( Lz .gt. 0. )taper = ( 1. - x ) * taper
            if( Lz .lt. -Lzcrt( icmp ) )taper = 0.
          end if
          done = .true.
c Miyamoto function - cannot reverse any Lz
        else if( cdft( icmp ) .eq. 'MIYA' )then
          done = .true.
c Sawamura's rule
        else if( cdft( icmp ) .eq. 'SAWA' )then
          if( dfcns( 2, icmp ) .ne. 0. )then
            x = Lz / Lzcrt( icmp )
            a = x
            ti = a
            do k = 1, 20
              a = a * ( 1. - x * x ) * real( 2 * k - 1 ) / real( 2 * k )
              ti = ti + a
            end do
            taper = taper * ( 1. + dfcns( 2, icmp ) * ti )
          end if
          done = .true.
        end if
      end if
c
      if( .not. done )call crash( 'TAPER',
     +                                 'Unknown distribution function' )
      return
      end
