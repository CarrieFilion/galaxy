      real*8 function vmean( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the expected mean orbital velocity for the given DF at the
c   requested radius
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
      real*8 Gammaf, vcirc, vmagris, vmlia
c
c local variables
      integer i, ifail, n
      real*8 a, arg, dfm, r1, r2
      include 'inc/pi.f'
c
      dfm = dfcns( 1, icmp )
      if( disc( icmp ) )then
c Kalnajs function
        if( ( cdft( icmp ) .eq. 'KALN' ) .and. kcold )then
          vmean = vcirc( r )
c Kalnajs function for Kuz'min/Toomre disc
        else if( ( cdft( icmp ) .eq. 'KALN' ) .and.
     +         ( ( ctype( icmp ) .eq. 'KT  ' ) .or.
     +           ( ctype( icmp ) .eq. 'SOFT' ) ) )then
          vmean = vmagris( r )
c generalised Kalnajs function
        else if( cdft( icmp ) .eq. 'LIA ' )then
          vmean = vmlia( r )
c Miyamoto's function for Kuz'min/Toomre disc
        else if( cdft( icmp ) .eq. 'MIYA' )then
          vmean = cmiy * r / ( 1. + r * r )**.75
c radius is a dummy argument - MTZ disc is self-similar
        else if( cdft( icmp ) .eq. 'ZANG' )then
          arg = .5 * ( dfm + 2. )
          ifail = 1
          vmean = sqrt( 2. / ( dfm + 1. ) ) * Gammaf( arg, ifail )
          if( ifail .eq. 0 )then
            arg = .5 * ( dfm + 1. )
            ifail = 1
            vmean = vmean / Gammaf( arg, ifail )
            if( ifail .eq. 0 )then
c dfm too large, mean v is very close to 1
              vmean = 1.
            end if
          end if
c Omega model
        else if( cdft( icmp ) .eq. 'OMEG' )then
c          vmean = sqrt( 3. * pi / 4. ) * omegak * r
          vmean = omegak * r
c Neil's function - use circular speed as a fudge
c        else if( cdft( icmp ) .eq. 'NEIL' )then
c          vmean = vcirc( r )
c Composite Omega model B0
        else if( cdft( icmp ) .eq. 'CPB0' )then
          vmean = r / sqrt( .75 * pi )
c Evans's DFs for power law discs - preprint March 1993
        else if( cdft( icmp ) .eq. 'EVAN' )then
          if( r .gt. 0. )then
            arg = 1 + .5 * dfm
            ifail = 0
            vmean = Gammaf( arg, ifail )
            arg = .5 + .5 * dfm
            vmean = vmean / Gammaf( arg, ifail )
            if( evbeta .lt. 0.d0 )then
              arg = -1.5 - ( 1 + dfm ) / evbeta
              vmean = vmean * Gammaf( arg, ifail )
              arg = -1 - ( 1 + dfm ) / evbeta
              vmean = vmean / Gammaf( arg, ifail )
            else
              arg = 2 + ( 1 + dfm ) / evbeta
              vmean = vmean * Gammaf( arg, ifail )
              arg = 2.5 + ( 1 + dfm ) / evbeta
              vmean = vmean / Gammaf( arg, ifail )
            end if
            vmean = vmean / sqrt( .5d0 * abs( evbeta ) * r**evbeta )
          else
            vmean = 0
          end if
c Evans/Collett DFs for Rybicki discs - MN v264 p353
        else if( cdft( icmp ) .eq. 'EVCO' )then
          n = nint( dfm )
          r2 = r * r
          r1 = sqrt( 1. + r2 )
c Appendix C
          if( n .ge. 0 )then
            vmean = 1
            a = 2**n * ( 1. + r2 ) / ( 1. + r1 )**n
            do i = 0, n
              if( i .gt. 0 )a =
     +           a * dble( 2 * i ) * r2 * dble( n + 1 - i )
     +                 / ( ( 1. + r1 ) * dble( 2 * i * ( 2 * i - 1 ) ) )
              vmean = vmean + a / sqrt( dble( n + i + 1 ) )
            end do
            vmean = sqrt( 2. / pi ) * vmean / ( r1 * ( 1. + r1 ) )
          else
            vmean = 0
          end if
        else
          call crash( 'VMEAN', 'Unknown type of DF' )
        end if
c non-rotating models
      else
c Merritt function for Jaffe halo
        if( cdft( icmp ) .eq. 'MERR' )then
          vmean = 0
c Dejonghe function for Plummer halo
        else if( cdft( icmp ) .eq. 'DEJO' )then
          vmean = 0
c Polyachenko-Shukhman function for uniform sphere
        else if( cdft ( icmp ) .eq. 'USPS' )then
          vmean = 0
        else
          call crash( 'VMEAN', 'Unknown type of DF' )
        end if
      end if
      return
      end
