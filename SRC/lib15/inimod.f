      subroutine inimod
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set variables in / model / after partial information has
c   been supplied from the input stream or header records
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/model.f'
c
c externals
      real*8 errfn, Gammaf
c
c local variables
      character*4 stype, type
      integer ifail, iuse, jcmp
      real w0
      real*8 arg, h, u
      save iuse
      include 'inc/pi.f'
      data iuse / 0 /
c
      jcmp = icmp
c assume false unless proved otherwise
      if( iuse .eq. 0 )then
        cmprssd = .false.
        fixed = .false.
        ldfit = .false.
        lgfld = .false.
        lhyb2 = .false.
        numcal = .false.
        singlr = .false.
        snglph = .false.
        iuse = 1
      end if
c value of iusepf flags whether to initiate a new dfiter run
      analyt = iusepf .eq. -1
      iusepf = 0
c work over all components
      do icmp = 1, ncmp
c set DF type flags
        dist( icmp ) = ( cdft( icmp ) .ne. 'NONE' ) .and.
     +                 ( cdft( icmp ) .ne. 'JEAN' )
        if( disc( icmp ) )then
          sphrod( icmp ) = .false.
c singular discs
          singlr = singlr .or. ( ctype( icmp ) .eq. 'MTZ ' ) .or.
     +                         ( ctype( icmp ) .eq. 'FMES' ) .or.
     +                         ( ( ctype( icmp ) .eq. 'POWE' ) .and.
     +                           ( dfcns( 3, icmp ) .lt. 0. ) )
c truncated discs
          trnctd( icmp ) = .not. ( ( ctype( icmp ) .eq. 'MFK ' ) .or.
     +                             ( ctype( icmp ) .eq. 'HUNT' ) .or.
     +                             ( ctype( icmp ) .eq. 'COSI' ) .or.
     +                             ( ctype( icmp ) .eq. 'FMES' ) .or.
     +                             ( ctype( icmp ) .eq. 'SPLY' ) .or.
     +                             ( ctype( icmp ) .eq. 'PPLY' ) .or.
     +                             ( ctype( icmp ) .eq. 'UNKN' ) )
          numcal = numcal .or. ( ctype( icmp ) .eq. 'ARBT' )
c special cases
          if( ( ctype( icmp ) .eq. 'SOFT' ) .or.
     +        ( ctype( icmp ) .eq. 'GAUS' ) )then
            if( epsiln .lt. 0.d0 )call crash( 'INIMOD',
     +                                     'nonsense value of epsilon' )
            if( ctype( icmp ) .eq. 'SOFT' )then
              if( epsiln .eq. 0.d0 )call crash( 'INIMOD',
     +                               'epsilon = 0 for a softened disk' )
              rstar = 1. + epsiln
            end if
c surface density normalization for power-law - see David Earn's thesis p 111
          else if( ctype( icmp ) .eq. 'POWE' )then
            if( abs( dfcns( 3, icmp ) ) .lt. 0.5 )then
              evbeta = -2. * dfcns( 3, icmp )
              arg = .5 * evbeta
              u = 1
              h = .5
              evsigma =
     +             Gammaf( u - arg, ifail ) * Gammaf( h + arg, ifail ) /
     +           ( Gammaf( u + arg, ifail ) * Gammaf( h - arg, ifail ) )
              evsigma = evsigma / ( 2. * pi )
            end if
          end if
        else
c finite halos
          trnctd( icmp ) = .not. ( ( type .eq. 'UNIS' ) .or.
     +                             ( type .eq. 'KING' ) .or.
     +                             ( type .eq. 'POLY' ) .or.
     +                             ( type .eq. 'UNKN' ) )
c singular halos
          singlr = singlr .or. ( ctype( icmp ) .eq. 'JAFF' ) .or.
     +                         ( ctype( icmp ) .eq. 'KEPL' ) .or.
     +                         ( ctype( icmp ) .eq. 'AGME' ) .or.
     +                         ( ctype( icmp ) .eq. 'HERN' ) .or.
     +                         ( ctype( icmp ) .eq. 'SISP' )
c fixed rotation curve models
          fixed = fixed .or. ( ctype( icmp ) .eq. 'FE  ' ) .or.
     +                       ( ctype( icmp ) .eq. 'BSS ' ) .or.
     +                       ( ctype( icmp ) .eq. '3198' ) .or.
     +                       ( ctype( icmp ) .eq. 'POWE' ) .or.
     +                       ( ctype( icmp ) .eq. 'POWC' )
c aspherical halos
          sphrod( icmp ) = ( ctype( icmp ) .eq. 'DFIT' ) .or.
     +                     ( ctype( icmp ) .eq. 'KKSP' ) .or.
     +                     ( ctype( icmp ) .eq. 'MYNG' )
c tabulated potential
          ldfit = ldfit .or. ctype( icmp ) .eq. 'DFIT'
          numcal = numcal .or. ldfit
          mnorm( icmp ) = 1
c setup isotropic models, preserving scaling
          if( ( ctype( icmp ) .eq. 'KING' ) .or.
     +        ( ctype( icmp ) .eq. 'POLY' ) )then
            call isoset
c normalization constant for Hernquist's cored Gaussian model, his eq (2.3)
          else if( ctype( icmp ) .eq. 'ISOG' )then
            arg = cmppar( 1, icmp )
c error function
            u = errfn( arg, ifail )
            u = 1.d0 - sqrt( pi ) * arg * exp( arg**2 ) * ( 1.d0 - u )
            cmppar( 2, icmp ) = 1.d0 / u
c normalization constant for Einasto halo mass
          else if( ctype( icmp ) .eq. 'EINA' )then
            h = dfcns( 3, icmp )
            arg = 3.d0 / h
            ifail = 0
            mnorm( icmp ) = ( .5d0 * h )**arg *
     +             exp( 2.d0 / h ) * Gammaf( arg, ifail ) / ( 4.d0 * h )
          else if( ldfit )then
c .pft file will be needed unless this is called from dfiter
            initls = .true.
            if( analyt )then
              if( impar( icmp ) .eq. 2 )then
                stype = ctype( icmp )
                w0 = dfcns( 3, icmp )
                ctype( icmp ) = 'KING'
                dfcns( 3, icmp ) = dfcns( 1, icmp )
                call isoset
                ctype( icmp ) = stype
                dfcns( 3, icmp ) = w0
              end if
            end if
          end if
        end if
      end do
c initialize distribution function constants
      do icmp = 1, ncmp
        if( dist( icmp ) )then
          call setcon
        end if
      end do
      icmp = jcmp
      return
      end
