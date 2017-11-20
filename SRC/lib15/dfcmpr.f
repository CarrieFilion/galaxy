      real*8 function dfcmpr( En, L )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      real*8 En, L
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 ajrtab, dfdejo, dfeddi, dfEina, dfhern, dfking, dfmerr
      real*8 dfNFW, E0tab
c
c local variables
      real*8 ajr, cnst, E0
c
      if( ctype( icmp ) .ne. 'ADIA' )call crash( 'DFCMPR',
     +                                               'Wrong halo type' )
c compute original energy for same radial action
      ajr = ajrtab( En, L )
      E0 = E0tab( ajr, L )
      dfcmpr = 0
      if( E0 .lt. Emaxo( icomp ) )then
        cnst = dfcns( 1, icmp )
        dfcns( 1, icmp ) = dfm0( icomp )
c control for uncompressed DF type
        if( idfn0( icomp ) .eq. 7 )then
c Merritt function for a Jaffe halo
          dfcmpr = dfmerr( E0, L )
        else if( idfn0( icomp ) .eq. 8 )then
c Dejonghe function for a Plummer halo
          dfcmpr = dfdejo( E0, L )
        else if( idfn0( icomp ) .eq. 12 )then
c King model
          dfcmpr = dfking( E0 )
        else if( idfn0( icomp ) .eq. 20 )then
c Hernquist function for a Hernquist halo
          dfcmpr = dfhern( E0, L )
        else if( idfn0( icomp ) .eq. 23 )then
c Eddington's formula for an isotropic spherical model
          if( ihalo0( icomp ) .eq. 25 )then
c pre-tablated version for an NFW halo
            dfcmpr = dfNFW( E0 )
          else if( ihalo0( icomp ) .eq. 32 )then
c pre-tablated version for an Einasto halo
            dfcmpr = dfEina( E0 )
          else
c generic version
            dfcmpr = dfeddi( E0 )
          end if
        else
          call crash( 'DFCMPR', 'Unrecognized orig DF type' )
        end if
        dfcns( 1, icmp ) = cnst
      end if
      return
      end
