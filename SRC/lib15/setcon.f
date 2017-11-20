      subroutine setcon
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c sets up values in / model / needed for computations of various DFs
c   called only from INIMOD
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
      include 'inc/setup.f'
c
c externals
      real*8 Gammaf, Gam1o2, Phitot, vcirc
c
c local valiables
      integer ifail, im, in, is, m2, m2s, ns
      logical check
      real a2n
      real*8 an, arg, a1, a2, a3, dfb, dfm
      include 'inc/pi.f'
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'SETCON',
     +                                        'Nonsense value of icmp' )
c set constants
      dfm = dfcns( 1, icmp )
      dfb = dfcns( 2, icmp )
c
      check = .false.
      arg = 0
      Phi0 = Phitot( arg )
      if( disc( icmp ) )then
        if( rscale( icmp ) .ne. 1. )call crash( 'SETCON',
     +                                            'rscale not equal 1' )
c Kalnajs's function
        if( cdft( icmp ) .eq. 'KALN' )then
          kcold = dfm .gt. 998.5
          if( ctype( icmp ) .eq. 'KT  ' )rstar = 1.
          if( ctype( icmp ) .eq. 'ISOC' )rstar = 2.
          if( ( rstar .le. 0.d0 ) .and.
     +        ( .not. kcold ) )call crash( 'SETCON',
     +                                  'r_* not set for a Kalnajs DF' )
          check = .true.
c Lia's function for Kuz'min/Toomre disc
        else if( cdft( icmp ) .eq. 'LIA ' )then
          check = ctype( icmp ) .eq. 'KT  '
c Miyamoto's function for Kuz'min/Toomre disc
        else if( cdft( icmp ) .eq. 'MIYA' )then
          check = ctype( icmp ) .eq. 'KT  '
          mmiy = dfm + .5
          nmiy = mmiy
          cmiy = ( real( 2 * mmiy ) - 1.5 ) / real( 2 * mmiy + 4 )
          cmiy = sqrt( cmiy )
          mmiy1 = mmiy + 1
          nmiy1 = nmiy + 1
c coefficients for even series
          amiy( 1, 1 ) = 3. / ( 16. * pi**2 )
          do im = 1, mmiy
            m2 = 2 * im
            amiy( 1, 1 ) = amiy( 1, 1 ) *
     +                     real( (m2 + 3) * (m2 + 2) * m2 ) /
     +                     real( 4 * im * (im + 2) * (m2 - 1) )
          end do
          do is = 2, mmiy1
            m2s = 2 * ( mmiy - is + 1 )
            amiy( 1, is ) = amiy( 1, is - 1 ) *
     +                    real( ( m2s + 1 ) * ( mmiy - is + 2 )**2 ) /
     +              real( ( is - 1 ) * ( m2s + 2 ) * ( mmiy + is + 1 ) )
          end do
c coefficients for odd series
          amiy( 2, 1 ) = 63. * cmiy / ( 2.**6.5 * pi**2 )
          do in = 1, nmiy
            an = in
            a2n = 2 * in
            amiy( 2, 1 ) = amiy( 2, 1 ) *
     +                    ( a2n + 4.5 ) * ( a2n + 3.5 ) * ( a2n + 2. ) /
     +                ( 4. * ( a2n + 1. ) * ( an + 1. ) * ( an + 2.5 ) )
          end do
          do is = 2, nmiy1
            ns = nmiy - is + 1
            amiy( 2, is ) = amiy( 2, is - 1 ) *
     +                real( ( 2 * ns + 3 ) * ( ns + 1 ) * ( ns + 2 ) ) /
     +                    ( real( ( 2 * ns + 4 ) * ( is - 1 ) ) *
     +                                     ( real( nmiy + is ) + 1.5 ) )
          end do
c DF for Mestel/Toomre/Zang disc
        else if( cdft( icmp ) .eq. 'ZANG' )then
          check = ctype( icmp ) .eq. 'MTZ '
c const = (m+1)**(m/2+1) / ( (m-1)/2 )!
          arg = .5 * dfm + .5
          an = .5 * dfm + 1.
          dfnorm( icmp ) = 1
          do while ( arg .gt. 2.d0 )
            arg = arg - 1
            dfnorm( icmp ) = dfnorm( icmp ) * ( dfm + 1. ) / arg
            an = an - 1
          end do
          ifail = 0
          dfnorm( icmp ) = dfnorm( icmp ) * ( dfm + 1. )**an /
     +                                              Gammaf( arg, ifail )
          ifail = 0
          arg = .5
          dfnorm( icmp ) = dfnorm( icmp ) / Gammaf( arg, ifail )
          dfnorm( icmp ) = dfnorm( icmp ) /
     +                                    ( 2.**( .5 * dfm + 1. ) * pi )
c Kalnajs Omega model
        else if( cdft( icmp ) .eq. 'OMEG'  )then
          check = ctype( icmp ) .eq. 'MFK '
          omegak = sqrt( 3. * pi / 4. ) * dfm
          omeg21 = 3. * pi / 4. - omegak * omegak
          dfnorm( icmp ) = 3. / ( 2. * pi * sqrt( omeg21 ) )
          Lzrmn( icmp ) = .5 * ( dfm - 1. ) * sqrt( 3. * pi / 4. )
c Kalnajs composite Omega model B0 - no parameter needed
        else if( cdft( icmp ) .eq. 'CPB0' )then
          check = ctype( icmp ) .eq. 'MFK '
          dfnorm( icmp ) = 2. / pi**2
          Lzrmn( icmp ) = -.5 * sqrt( 3. * pi / 4. )
c Evans-Collett functions for the Rybicki disc (MNRAS v264, p353)
        else if( cdft( icmp ) .eq. 'EVCO' )then
          check = ctype( icmp ) .eq. 'RYBI'
          kcold = dfm .gt. 998.5
c Sawamura's functions for polynomial discs (PASJ v40, p279)
        else if( cdft( icmp ) .eq. 'SAWA' )then
          check = ctype( icmp ) .eq. 'SPLY'
          Lzcrt( icmp ) = vcirc( 1.d0 )
c disc function by Lucy method - parameters already defined in dfm & dfb
c        else if( cdft( icmp ) .eq. 'LUCY' )then
c Neil's disc set-up
c        else if( cdft( icmp ) .eq. 'NEIL' )then
c          dfcns( 1, icmp ) = dfm
c          rtrunc( icmp ) = dfb
c Evans functions for power law discs (Evans & Read 1998, MNRAS, v300, p83)
        else if( cdft( icmp ) .eq. 'EVAN' )then
          check = ctype( icmp ) .eq. 'POWE'
          if( abs( evbeta ) .gt. 1.d0 )call crash( 'SETCON',
     +                       'Impossible power law index for EVAN fn' )
          dfb = ( 1 + dfm ) / evbeta
          dfcns( 2, icmp ) = dfb
c declining rotation curves
          if( evbeta .gt. 0.d0 )then
            if( dfm .le. -1.d0 )call crash( 'SETCON',
     +                              'Impossible DF index for EVAN fn' )
c \Gamma(2+(1+q)/\beta) / ( \Gamma((1+q)/2) * \Gamma(1+(1+q)/\beta-q/2) )
            a1 = 2 + ( 1 + dfm ) / evbeta
            a2 = .5 * ( 1 + dfm )
            a3 = 1 + ( 1 + dfm ) / evbeta - .5 * dfm
            dfnorm( icmp ) = Gam1o2( a1, a2, a3 )
          else
c rising rotation curves
            if( dfm .le. -evbeta - 1.d0 )call crash( 'SETCON',
     +                               'Impossible DF index for EVAN fn' )
c \Gamma(q/2-(1+q)/\beta) / ( \Gamma((1+q)/2) * \Gamma(-1-(1+q)/\beta) )
            a1 = .5 * dfm - ( 1 + dfm ) / evbeta
            a2 = .5 * ( 1 + dfm )
            a3 = -1 - ( 1 + dfm ) / evbeta
            dfnorm( icmp ) = Gam1o2( a1, a2, a3 )
          end if
c additional normalising factors
          arg = .5 * dfm
          dfnorm( icmp ) = dfnorm( icmp ) * evsigma /
     +                                        ( 2.d0**arg * sqrt( pi ) )
c I choose units such that \psi_a = 1 / \beta, i.e. v_\beta = 1
          dfnorm( icmp ) = dfnorm( icmp ) * abs( evbeta )
c the following factor is now included in the DF
c          arg = ( 1 + dfm ) / evbeta
c          dfnorm( icmp ) = dfnorm( icmp ) * abs( evbeta )**arg
c Polyachenko's lowered polytropes for discs (Astr lett 1994 v20, p416)
        else if( cdft( icmp ) .eq. 'PPLY' )then
          check = ctype( icmp ) .eq. 'PPLY'
          arg = .5 * ( dfb + 1. )
          dfnorm( icmp ) = 2.**( .5 * dfb ) * Gammaf( arg, ifail )
          arg = .5
          dfnorm( icmp ) = dfnorm( icmp ) * Gammaf( arg, ifail )
          arg = dfm + 1.
          dfnorm( icmp ) = dfnorm( icmp ) * Gammaf( arg, ifail )
          arg = dfm + .5 * dfb + 2.
          dfnorm( icmp ) = dfnorm( icmp ) / Gammaf( arg, ifail )
        end if
      else
c Osipkov-Merritt function for a Jaffe halo
        if( cdft( icmp ) .eq. 'MERR' )then
          check = ctype( icmp ) .eq. 'JAFF'
c anisotropy radius is infinite for isotropic model, flagged as dfm -ve
          dfb = 0.
          if( dfm .gt. 0. )dfb = 1. / dfm**2
c Dejonghe function for a Plummer halo
        else if( cdft( icmp ) .eq. 'DEJO' )then
          check = ctype( icmp ) .eq. 'PLUM'
c q =< 2
          if( dfm .gt. 2. )call crash( 'SETCON',
     +                              'Dejonghe q outside allowed range' )
c King function
        else if( cdft( icmp ) .eq. 'KING' )then
          check = ctype( icmp ) .eq. 'KING'
c spherical polytrope
        else if( cdft( icmp ) .eq. 'POLY' )then
          check = ctype( icmp ) .eq. 'POLY'
c Kuz'min-Kutuzov functions for oblate spheroid
        else if( cdft( icmp ) .eq. 'DJDZ' )then
          check = ctype( icmp ) .eq. 'KKSP'
c Hernquist's function
        else if( cdft( icmp ) .eq. 'HERN' )then
          check = ctype( icmp ) .eq. 'HERN'
c Henon's isotropic function for the spherical isochrone
        else if( cdft( icmp ) .eq. 'HNON' )then
          check = ctype( icmp ) .eq. 'ISOC'
c Eddington's formula for an isotropic DF for a spherical system
        else if( cdft( icmp ) .eq. 'EDDI' )then
          check = .not. sphrod( icmp )
c singular isothermal sphere
        else if( cdft( icmp ) .eq. 'SISP' )then
          check = ctype( icmp ) .eq. 'SISP'
c Polyachenko-Shukhman DF for uniform sphere - see BT08 problem 4.11
        else if( cdft( icmp ) .eq. 'USPS' )then
          check = ctype( icmp ) .eq. 'UNIS'
          dfnorm( icmp ) = 0.75 / ( pi**3 * rscale( icmp )**2 )
        end if
      end if
c none of the above are for composite models
      check = check .and. ( ncmp .eq. 1 )
c numerically determined potential from DFITER
      if( cdft( icmp ) .eq. 'DFIT' )then
        check = ctype( icmp ) .eq. 'DFIT'
        if( impar( icmp ) .eq. 2 )sigmai2 =
     +                       ( dfcns( 3, icmp ) - Phitot( 0.d0 ) ) / dfm
c numerically determined energies from adiabatic compression
      else if( cdft( icmp ) .eq. 'COMP' )then
        check = ctype( icmp ) .eq. 'ADIA'
        if( impar( icmp ) .eq. 2 )sigmai2 =
     +                       ( dfcns( 3, icmp ) - Phitot( 0.d0 ) ) / dfm
c Shu's epicyclic DF for generic disks
      else if( cdft( icmp ) .eq. 'SHUE' )then
        check = .true.
      end if
c warn if potential and DF are inconsistent
      if( .not. check )then
        if( master )print
     +            '( '' Disc/halo type = '', a4, '' DF type = '', a4 )',
     +                                       ctype( icmp ), cdft( icmp )
        if( ncmp .eq. 1 )then
          call crash( 'SETCON', 'Unrecognized DF for this potential' )
        else if( master )then
          print *, 'Warning: potential and DF not consistent'
          write( no, * )'Warning: potential and DF not consistent'
        end if
      end if
      return
      end
