      real*8 function sglint( Lz1, Lz2, E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c integrates a two-integral DF for an axisymmetric model or disk over the
c   given range of angular momenta and at the constant given energy
c may use QUADPACK routine DQAGS
c
c calling arguments
      real*8 E, Lz1, Lz2
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      integer mc, nc
      parameter ( nc = 11, mc = 101 )
      real*8 Ecut, Lzc( mc ), mpart( 2 * nc ), ord( mc )
      common / cut / Ecut, Lzc, mpart, ord
c
      include 'inc/model.f'
c
c externals
      external dfnje
      real*8 dfnje, Lzmax, Lzmin, quad_tab, quad_gnr
c
c local arrays
      real*8 absc( mc ), val( mc )
c
c local variables
      integer ic, ier, ip, iuse, jc, kc
      logical calc, lbound, ubound
      real*8 epsa, epsr, Er, fract, Lz, Lzm
      save iuse
      include 'inc/pi.f'
c
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
        Ecut = 1.d+10
        iuse = 1
      end if
c integral analytic for Omega model
      if( cdft( icmp ) .eq. 'OMEG' )then
        fract = max( omeg21 - 2. * ( E - Phi0 - omegak * Lz2 ), 0.d0 )
        sglint = sqrt( fract )
        fract = max( omeg21 - 2. * ( E - Phi0 - omegak * Lz1 ), 0.d0 )
        sglint = sglint - sqrt( fract )
        sglint = sqrt( 4. * pi / 3. ) * dfnorm( icmp ) * sglint / omegak
c also for Polyachenko-Shukhman DF for the uniform sphere
c   see herschel:/home/sellwood/docs/models/harmonic.tex
      else if( cdft( icmp ) .eq. 'USPS' )then
        fract = max( Lz2 * Lz2 - Phimax - 2. * ( E - Phi0 ), 0d0 )
        sglint = sqrt( fract )
        fract = max( Lz1 * Lz1 - Phimax - 2. * ( E - Phi0 ), 0d0 )
        sglint = sglint - sqrt( fract )
        sglint = 8. * pi**3 * dfnorm( icmp ) *
     +               sqrt( rscale( icmp )**3 / cmpmas( icmp ) ) * sglint
c spheroids
      else if( sphrod( icmp ) )then
c avoid problems caused by an effectively zero range of integration
        if( Lz2 - Lz1 .lt. 1.d-8 )then
          sglint = 0
        else
c clear table if value of energy has changed
          if( E .ne. Ecut )then
            Lz = Lzmin( E )
            Lzm = Lzmax( E )
            Ecut = E
            do ic = 1, mc
              if( ic .eq. 1 )then
                Lzc( 1 ) = Lz
              else if( ic .eq. mc )then
                Lzc( mc ) = Lzm
              else
                Lzc( ic ) = Lz +
     +                    real( ic - 1 ) * ( Lzm - Lz ) / real( mc - 1 )
              end if
              ord( ic ) = dfnje( Lzc( ic ) )
            end do
          end if
c create table of values for numerical integration
          ic = 1
          jc = 0
c from possibly larger than Lzmin
          if( Lz1 .gt. Lzc( 1 ) )then
            jc = 1
            absc( 1 ) = Lz1
            val( 1 ) = dfnje( Lz1 )
            do while ( Lz1 .lt. Lzc( ic ) )
              ic = ic + 1
              if( ic .gt. mc )call crash( 'SGLINT', 'Logical error' )
            end do
          end if
c copy usable values
          do kc = ic, mc
            if( Lz2 .ge. Lzc( kc ) )then
              jc = jc + 1
              absc( jc ) = Lzc( kc )
              val( jc ) = ord( kc )
            end if
          end do
c to possibly less than Lzmax
          kc = min( ic + jc + 1, mc )
          if( Lz2 .lt. Lzc( kc ) )then
c an extra abscissa only when sufficiently far beyond the last
            if( abs( absc( jc ) - Lz2 ) .gt. .01 * Lz2 )jc = jc + 1
            absc( jc ) = Lz2
            val( jc ) = dfnje( Lz2 )
          end if
c evaluate integral
          if( jc .gt. 3 )then
            sglint = quad_tab( absc, val, jc, er )
          else
c trapezium rule will do
            sglint = 0
            do ic = 2, jc
              sglint = sglint + .5 * ( val( ic ) + val( ic -1 ) ) *
     +                              ( absc( ic ) - absc( ic -1 ) )
            end do
          end if
c double value to include retrograde part
          sglint = 2 * sglint
        end if
      else
c clear table if value of energy has changed
        if( E .ne. Ecut )then
          Ecut = E
          Lz = Lzmin( E )
          Lzm = Lzmax( E )
          do ic = 1, nc
            Lzc( ic + nc ) = Lz
     +                  + real( ic - 1 ) * ( Lzm - Lz ) / real( nc - 1 )
            mpart( ic + nc ) = -1
c retrograde range
            Lzc( nc - ic + 1 ) = -Lzc( ic + nc )
            mpart( nc - ic + 1 ) = -1
          end do
        end if
c sum over retrograde then direct parts
        sglint = 0
        do ip = 1, 2
          do ic = 1, nc - 1
            jc = nc * ( ip - 1 ) + ic
c skip if part is totally outside Lz range requested
            if( .not. ( ( Lz1 .gt. Lzc( jc + 1 ) ) .or.
     +                  ( Lz2 .lt. Lzc( jc )     ) ) )then
              lbound = Lz1 .le. Lzc( jc )
              ubound = Lz2 .ge. Lzc( jc + 1 )
c determine whether integration is required
              calc = .not. ( lbound .and. ubound )
              calc = calc .or. ( mpart( jc ) .lt. 0. )
              if( calc )then
                Lz = Lz1
                if( lbound )Lz = Lzc( jc )
                Lzm = Lz2
                if( ubound )Lzm = Lzc( jc + 1 )
c avoid problems caused by an effectively zero range of integration
                if( Lzm - Lz .lt. 1.d-8 )then
                  fract = 0
                else
                  epsa = 1.d-8
                  epsr = epsa
c adaptive integrator
                  fract = quad_gnr( dfnje, Lz, Lzm, epsa, epsr, ier )
                  if( ier .ne. 0 )then
                    epsa = abs( ( Lzm - Lz ) / Lzm )
                    print *, 'ier =', ier,
     +          ' from QUAD_GNR for fraction', sngl( epsa ), ' of range'
c                    if(
c     +         epsa .gt. 1.d-4 )call crash( 'SGLINT', 'QUADPACK error' )
                  end if
                end if
                sglint = sglint + fract
c save result if it is re-usable
                if( lbound .and. ubound )mpart( jc ) = fract
              else
                sglint = sglint + mpart( jc )
              end if
            end if
c end loop over fractions
          end do
c end loop over retrograde & direct parts
        end do
      end if
      return
      end
