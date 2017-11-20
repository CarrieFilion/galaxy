      real*8 function distfn( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c driver routine for all programmed DF forms
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      real*8 dfcmpr, dfcpb0, dfDejo, dfDJDZ, dfEddi, dfEina, dfEvan
      real*8 dfEvCo, dfHern, dfHnon, dfitfn, dfKaln, dflia,  dfKing
      real*8 dfMerr, dfMiym, dfMTZ,  dfNFW,  dfOmeg, dfPoly, dfSawa
      real*8 dfShu,  dfunis
      real*8 dfwght, taper
c, dflucy
c
c local variables
      logical ltaper
      real*8 sfac
      include 'inc/pi.f'
c
      ltaper = .false.
      if( disc( icmp ) )then
        sfac = cmpmas( icmp )
c Kalnajs function
        if( cdft( icmp ) .eq. 'KALN' )then
          distfn = dfkaln( E, Lz )
          ltaper = .true.
c Lia's function
        else if( cdft( icmp ) .eq. 'LIA ' )then
          distfn = dflia( E, Lz )
          ltaper = .true.
c Miyamoto function
        else if( cdft( icmp ) .eq. 'MIYA' )then
          distfn = dfmiym( E, Lz )
          ltaper = .true.
c Zang function
        else if( cdft( icmp ) .eq. 'ZANG' )then
          if( Lz .eq. 0. )then
            distfn = 0
          else
            distfn = dfmtz( E, Lz )
          end if
          ltaper = .true.
c Kalnajs Omega model
        else if( cdft( icmp ) .eq. 'OMEG' )then
          distfn = dfomeg( E, Lz )
c Lucy type DF for a disc
c        else if( cdft( icmp ) .eq. 'LUCY' )then
c          distfn = dflucy( E, Lz )
c Composite Omega model B0
        else if( cdft( icmp ) .eq. 'CPB0' )then
          distfn = dfcpb0( E, Lz )
c Wyn Evans's DF for power law discs
        else if( cdft( icmp ) .eq. 'EVAN' )then
          distfn = dfevan( E, Lz )
          ltaper = .true.
c Evans-Collett DF for the Rybicki disc
        else if( cdft( icmp ) .eq. 'EVCO' )then
          distfn = dfevco( E, Lz )
          ltaper = .true.
c Sawamura's DF for his polynimial disc
        else if( cdft( icmp ) .eq. 'SAWA' )then
          distfn = dfsawa( E, Lz )
          ltaper = .true.
c Shu's epicyclic DF for a disc
        else if( cdft( icmp ) .eq. 'SHUE' )then
          distfn = dfShu( E, Lz )
          ltaper = .false.
        else
          call crash( 'DISTFN', 'Unknown distribution function type' )
        end if
      else
        sfac = 1
c Merritt function for a Jaffe halo
        if( cdft( icmp ) .eq. 'MERR' )then
          distfn = dfmerr( E, Lz )
c Dejonghe function for Plummer halo
        else if( cdft( icmp ) .eq. 'DEJO' )then
          distfn = dfdejo( E, Lz )
c King model
        else if( cdft( icmp ) .eq. 'KING' )then
          distfn = dfking( E )
c Uniform Sphere
        else if( cdft( icmp ) .eq. 'USPS' )then
          distfn = dfunis( E, Lz )
c spherical polytrope
        else if( cdft( icmp ) .eq. 'POLY' )then
          distfn = dfpoly( E )
c Dejonghe-deZeeuw function for Kuz'min-Kutuzov spheroid
        else if( cdft( icmp ) .eq. 'DJDZ' )then
          distfn = dfdjdz( E, Lz )
          ltaper = .true.
c isotropic form for DF iteration
        else if( cdft( icmp ) .eq. 'DFIT' )then
          distfn = dfitfn( E )
c Hernquist function for a Hernquist halo
        else if( cdft( icmp ) .eq. 'HERN' )then
          distfn = dfhern( E, Lz )
c Henon's isotropic function for a spherical isochrone
        else if( cdft( icmp ) .eq. 'HNON' )then
          distfn = dfhnon( E )
c Eddington's formula for an isotropic spherical model
        else if( cdft( icmp ) .eq. 'EDDI' )then
          if( ctype( icmp ) .eq. 'NFW ' )then
c special version for an NFW halo
            distfn = dfNFW( E )
          else if( ctype( icmp ) .eq. 'EINA' )then
c special version for an Einasto halo
            distfn = dfEina( E )
          else
c generic version
            distfn = dfeddi( E )
          end if
c singular isothernal sphere - mass and radial scales /= 0 make no sense
        else if( cdft( icmp ) .eq. 'SISP' )then
          distfn = exp( -2. * E ) / ( 4. * pi**2.5 )
c adiabatically compressed halo
        else if( cdft( icmp ) .eq. 'COMP' )then
          distfn = dfcmpr( E, Lz )
        else
          call crash( 'DISTFN', 'Unknown distribution function type' )
        end if
      end if
c scale by mass of component for most disk types
      distfn = sfac * distfn
c weight distfn if required
      if( uqmass )distfn = distfn * dfwght( E, Lz )
c add angular momentum taper and retrograde star rule
      if( ltaper )distfn = distfn * taper( Lz )
      return
      end
