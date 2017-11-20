      subroutine msetup
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine for selection of disc and/or halo and DF types and other
c   parameters of a model
c values are read either from stdin in an interactive manner, or
c   from the runXXX.dat file if common variable ni is greater than zero
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      character dftype*4, dtype*4, htype*4
      integer lnblnk
      logical gtlogl
      real*8 gmassd, gmassh, gmassi, vcirc
c
c local arrays
      integer mdfn, mdiscs
      parameter ( mdfn = 28, mdiscs = 22 )
      integer npars( mdfn )
c
c local variables
      character line*120, type*4
      integer i, j, jcomp, m, n
      logical froz, local, lskip, term
      real cc
      real*8 r2, rtap
c
c number of parameters needed for the various DF types
      data npars / 0, 1, 2, 1, 1, 1, 1, 1, 2, 2,
     +             0, 0, 0, 1, 0, 1, 0, 1, 2, 1,
     +             1, 0, 0, 0, 0, 1, 0, 0 /
c
c determine type of input stream
      term = ni .le. 0
c check whether already set
      if( ncmp .gt. 0 )then
        if( master )print
     +                  *, 'model already set with', ncmp, ' components'
        if( gtlogl( 'Do you want to abort' ) )call crash( 'MSETUP',
     +                                 'Did not wish to re-initialize' )
      end if
c initialize
      type = dtype( 1 )
      if( ndiscs .gt. mdiscs )call crash( 'MSETUP', 'New disk type' )
      type = htype( 1 )
      if( nhalos .gt. 34 )call crash( 'MSETUP', 'New halo type' )
      type = dftype( 1 )
      if( ndfns .gt. mdfn )call crash( 'MSETUP', 'New DF type' )
      cmprssd = .false.
      jcomp = -1
c choose number of components
      if( term )then
        call gtintg( 'Enter number of mass components', ncmp )
        if( master )print *, 'Setting up', ncmp, ' components'
      else
        call getline( ni, line )
        if( line( 1:4 ) .ne. 'NCOM' )go to 1
        n = lnblnk( line )
        read( line( 11:n ), *, err = 1 )ncmp
        if( master )write( no, * )'Setting up', ncmp, ' components'
      end if
      if( ( ncmp .lt. 1 ) .or. ( ncmp .gt. mcmp ) )call crash( 'MSETUP',
     +                               'Impossible number of components' )
      froz = .false.
c work over components
      do icmp = 1, ncmp
        if( term )then
          if( master )print *, 'Starting component', icmp
          disc( icmp ) = gtlogl( 'Is this a disk' )
        else
          if( master )then
            write( no, * )
            write( no, * )'Starting component', icmp
          end if
          call getline( ni, line )
          if( line( 1:4 ) .ne. 'DISC' )go to 1
          n = lnblnk( line )
          read( line( 11:n ), *, err = 1 )disc( icmp )
c read disc/halo type
          call getline( ni, line )
          if( line( 1:4 ) .ne. 'TYPE' )go to 1
          n = lnblnk( line )
          m = 11
          do while ( line( m:m ) .eq. ' ' )
            m = m + 1
            if( m .gt. n )go to 1
          end do
          type = line( m:m+3 )
        end if
c disk
        if( disc( icmp ) )then
          imtyp( icmp ) = 0
          do while ( imtyp( icmp ) .eq. 0 )
            if( term )then
              call gtchar( 'Enter disc type', type )
              call uppercase( type )
            end if
            if( type( 1:3 ) .eq. 'END' )go to 2
c check for recognised name
            do i = 2, ndiscs
              if( type .eq. dtype( i ) )imtyp( icmp ) = i
            end do
            if( imtyp( icmp ) .eq. 0 )then
              if( master )print *, 'Unrecognized disc type'
              if( .not. term )go to 1
            end if
          end do
          ctype( icmp ) = type
          if( ( .not. term ) .and. master )write( no,
     +             '( '' Disc type selected is '', a4, '' index'', i3 )'
     +                                     )ctype( icmp ), imtyp( icmp )
c pick out finite disks
          trnctd( icmp ) = .not. ( ( type .eq. 'MFK ' ) .or.
     +                             ( type .eq. 'HUNT' ) .or.
     +                             ( type .eq. 'COSI' ) .or.
     +                             ( type .eq. 'FMES' ) .or.
     +                             ( type .eq. 'SPLY' ) .or.
     +                             ( type .eq. 'PPLY' ) )
          if( .not. trnctd( icmp ) )rtrunc( icmp ) = 1
c set defaults
          epsiln = 0
          impar( icmp ) = 0
          rhole = 0
          cc = 0
          rstar = 1
          if( ( .not. term ) .and.
     +        ( ( type .eq. 'SOFT' ) .or. ( type .eq. 'GAUS' ) .or.
     +          ( type .eq. 'POWE' ) .or. ( type .eq. 'COSI' ) .or.
     +          ( type .eq. 'SPLY' ) .or. ( type .eq. 'HUNT' ) ) )then
            m = m + 4
            if( type .eq. 'HUNT' )then
              read( line( m:n ), *, err = 1 )impar( icmp )
            else
              read( line( m:n ), *, err = 1 )cc
            end if
          end if
          if( ( type .eq. 'SOFT' ) .or. ( type .eq. 'GAUS' ) )then
            if( term )then
              call gtreal( 'Enter softening length', cc )
            else
              if( master )write( no,
     +                          '( '' Softening parameter'', f7.4 )' )cc
            end if
            epsiln = cc
            dfcns( 3, icmp ) = cc
          else if( type .eq. 'POWE' )then
            if( term )then
              call gtreal( 'Enter power law index', cc )
            else
              if( master )write( no,
     +                              '( '' Power law index'', f7.4 )' )cc
            end if
            dfcns( 3, icmp ) = cc
            if( cc .lt. -1.99 )then
              rhole = 1
            else
              rhole = .01
            end if
            if( abs( cc ) .lt. 1.d-3 )call crash( 'MSETUP',
     +                               'Unsuitable index, use dtype MTZ' )
          else if( type .eq. 'COSI' )then
            if( term )then
              call gtreal( 'Enter aspect ratio', cc )
            else
              if( master )write( no,
     +                    '( '' Aspect ratio of COSI disc'', f7.4 )' )cc
            end if
            if( ( cc .le. 0. ) .or.
     +          ( cc .gt. 1. ) )call crash( 'MSETUP',
     +                      'Impossible aspect ratio for cosine disc' )
            rhole = 1 - cc
            rtrunc( icmp ) = 1 + cc
            dfcns( 3, icmp ) = 1. / cc
          else if( type .eq. 'SPLY' )then
            if( term )then
              call gtreal( 'Enter index for Sawamura polynomial', cc )
            else
              if( master )write( no,
     +                      '( '' Index for Sawamura disc'', f7.4 )' )cc
            end if
            dfcns( 3, icmp ) = cc
          else if( type .eq. 'HUNT' )then
            if( term )then
              call gtintg(
     +                  'Enter index of Hunter disc', impar( icmp ) )
            else
              if( master )write( no,
     +               '( '' Index for Hunter disc'', i4 )' )impar( icmp )
            end if
            if( ( impar( icmp ) .lt. 1 ) .or.
     +          ( impar( icmp ) .gt. 5 ) )call crash( 'MSETUP',
     +                                'Invalid index for Hunter disc' )
          else if( ( type .eq. 'MTZ ' ) .or. ( type .eq. 'FMES' ) )then
            rhole = .01
          end if
        else
c halo
          imtyp( icmp ) = 0
          do while ( imtyp( icmp ) .eq. 0 )
            if( term )then
              call gtchar( 'Enter halo type', type )
              call uppercase( type )
            end if
            if( type( 1:3 ) .eq. 'END' )go to 2
c check for recognised name
            do i = 2, nhalos
              if( type .eq. htype( i ) )imtyp( icmp ) = i
            end do
            if( imtyp( icmp ) .eq. 0 )then
              if( master )print *, 'Unrecognized halo type'
              if( .not. term )go to 1
            end if
          end do
          if( type .eq. 'ASDI' )froz = .true.
          ctype( icmp ) = type
          if( ( .not. term ) .and. master )write( no,
     +             '( '' Halo type selected is '', a4, '' index'', i3 )'
     +                                     )ctype( icmp ), imtyp( icmp )
c finite halos
          trnctd( icmp ) = .not. ( ( type .eq. 'UNIS' ) .or.
     +                             ( type .eq. 'KING' ) .or.
     +                             ( type .eq. 'POLY' ) .or.
     +                             ( type .eq. 'UNKN' ) .or.
     +                             ( type .eq. 'CUBI' ) )
          if( .not. trnctd( icmp ) )rtrunc( icmp ) = 1
c fixed rotation curve models
          fixed = fixed .or. ( ctype( icmp ) .eq. 'FE  ' ) .or.
     +                       ( ctype( icmp ) .eq. 'BSS ' ) .or.
     +                       ( ctype( icmp ) .eq. '3198' ) .or.
     +                       ( ctype( icmp ) .eq. 'POWE' ) .or.
     +                       ( ctype( icmp ) .eq. 'POWC' )
c get halo constant
          if( ( .not. term ) .and.
     +        ( ( type .eq. 'FE  ' ) .or. ( type .eq. 'POWE' ) .or.
     +          ( type .eq. 'KING' ) .or. ( type .eq. 'POLY' ) .or.
     +          ( type .eq. 'ISOT' ) .or. ( type .eq. 'KKSP' ) .or.
     +          ( type .eq. 'AGME' ) .or. ( type .eq. 'DFIT' ) .or.
     +          ( type .eq. 'MYNG' ) .or. ( type .eq. 'ISOG' ) .or.
     +          ( type .eq. 'EINA' ) ) )then
            m = m + 4
            if( type .eq. 'DFIT' )then
              read( line( m:n ), *, err = 1 )impar( icmp )
            else
              read( line( m:n ), *, err = 1 )cc
            end if
          end if
          if( ( type .eq. 'FE  ' ) .or. ( type .eq. 'POWE' ) .or.
     +        ( type .eq. 'KING' ) .or. ( type .eq. 'POLY' ) .or.
     +        ( type .eq. 'ISOT' ) .or. ( type .eq. 'KKSP' ) .or.
     +        ( type .eq. 'AGME' ) .or. ( type .eq. 'MYNG' ) .or.
     +        ( type .eq. 'ISOG' ) .or. ( type .eq. 'EINA' ) )then
            if( term )then
              call gtreal( 'Enter halo const', cc )
            else
              if( master )write( no, '( '' Halo constant'', f7.4 )' )cc
            end if
            if( type .eq. 'ISOT' )cc = cc**2
          else if( type .eq. 'DFIT' )then
c extra parameter for DFIT type halos
            if( term )then
              call gtintg(
     +                'Enter 1 for polytrope, 2 for lowered isothermal',
     +                                                   impar( icmp ) )
            else
              if( master )write( no,
     +                           '( '' DFIT type'', i4 )' )impar( icmp )
            end if
c iusepf may be already set, and is preset to -1 for cold starts in DFITER only
            if( abs( iusepf ) .ne. 1 )iusepf = 0
          else if( ( type .eq. 'ADIA' ) .or. ( type .eq. 'VKA1' ) )then
            if( type .eq. 'ADIA' )imtyp( icmp ) = 23
            if( type .eq. 'VKA1' )imtyp( icmp ) = 27
            jcomp = icmp
          else
            cc = 0
          end if
          dfcns( 3, icmp ) = cc
          cmppar( 1, icmp ) = cc
          if( ( type .eq. 'KEPL' ) .and. ( rhole .eq. 0.d0 ) )rhole = 1
c check for sensible power laws
          if( type .eq. 'POWE' )then
            if( cc .lt. 0.0 )rhole = 1
            if( cc .lt. -0.5 )call crash( 'MSETUP',
     +                 'Unphysical index for power law rotation curve' )
c King model or polytrope - work in natural units
          else if( type .eq. 'KING' )then
            if( cc .le. 0. )call crash( 'MSETUP',
     +                              'Nonsense W0 value for King model' )
          else if( type .eq. 'POLY' )then
            if( ( cc .le. 0.5 ) .or. ( cc .ge. 5.0 )
     +         )call crash( 'MSETUP', 'Impossible index for polytrope' )
          else if( type .eq. 'AGME' )then
            if( ( cc .ge. 0. ) .or.
     +          ( cc .le. -3. ) )call crash( 'MSETUP',
     +                   'Unreasonable index for Aguilar-Merritt halo' )
c check for allowed range of polytropes
          else if( type .eq. 'POLY' )then
            if( cc .lt. 0.5d0 )call crash( 'MSETUP',
     +                                'Impossible index for polytrope' )
            if( cc .eq. 5.0d0 )call crash( 'MSETUP',
     +                  'Polytrope index 5 is same as a Plummer model' )
            if( cc .gt. 5.0d0 )call crash( 'MSETUP',
     +                                'Index too large for polytrope ' )
c reset constants for NGC 3198 halo
          else if( type .eq. '3198' )then
            call crash( 'MSETUP', '3198 halo no longer available' )
c            haloc = .63
c            hrad = 4.7
c reset constants for the TEDS halo
          else if( type .eq. 'TEDS' )then
            call crash( 'MSETUP', 'TEDS halo no longer available' )
c            haloc = 1.8 / 24.
c            hrad = 24. / 6.9
          end if
        end if
c mass
        if( type .ne. 'ISOT' )then
          if( term )then
           call gtreal( 'Enter mass of this component', cmpmas( icmp ) )
          else
            call getline( ni, line )
            if( line( 1:4 ) .ne. 'MASS' )go to 1
            n = lnblnk( line )
            read( line( 11:n ), *, err = 1 )cmpmas( icmp )
            if( master )write( no,
     +           '( '' Mass of this component'', f7.4 )' )cmpmas( icmp )
          end if
          if( cmpmas( icmp ) .lt. 0. )call crash( 'MSETUP',
     +                                'Nonsense value of mass entered' )
        end if
c radial scale
        if( term )then
          call gtreal(
     +          'Enter radial scale of this component', rscale( icmp ) )
        else
          call getline( ni, line )
          if( line( 1:4 ) .ne. 'SCAL' )go to 1
          n = lnblnk( line )
          read( line( 11:n ), *, err = 1 )rscale( icmp )
          if( master )write( no,
     +   '( '' Radial scale of this component'', f7.2 )' )rscale( icmp )
        end if
        if( rscale( icmp ) .le. 0. )call crash( 'MSETUP',
     +                        'Nonsense value of radial scale entered' )
        if( type .eq. 'ISOT' )cmpmas( icmp ) = cc * rscale( icmp )
c truncation radius
        if( trnctd( icmp ) )then
          if( term )then
            call gtreal(
     +      'Enter truncation radius in scale lengths', rtrunc( icmp ) )
          else
            read( line( 11:n ), *, err = 1 )cc, rtrunc( icmp )
            if( master )write( no,
     +   '( '' Truncation radius'', f7.3, '' scales'' )' )rtrunc( icmp )
          end if
          trnctd( icmp ) = ( rtrunc( icmp ) .lt. 999. ) .and.
     +                     ( type .ne. 'DFIT' )
        end if
        rtrunc( icmp ) = rscale( icmp ) * rtrunc( icmp )
c DF type
        if( .not. term )then
          call getline( ni, line )
          if( line( 1:4 ) .ne. 'DFTY' )go to 1
          n = lnblnk( line )
          m = 11
          do while ( line( m:m ) .eq. ' ' )
            m = m + 1
            if( m .gt. n )go to 1
          end do
          if( m .gt. n - 3 )go to 1
          type = line( m:m + 3 )
        end if
        idftyp( icmp ) = 0
        do while ( idftyp( icmp ) .eq. 0 )
          if( term )then
            call gtchar( 'Enter distribution function type', type )
            call uppercase( type )
          end if
          if( type( 1:3 ) .eq. 'END' )go to 2
c check for recognised name
          do i = 1, ndfns
            if( type .eq. dftype( i ) )idftyp( icmp ) = i
          end do
          if( idftyp( icmp ) .eq. 0 )then
            if( master )print *, 'Unrecognized DF type'
            if( .not. term )go to 1
          end if
        end do
        cdft( icmp ) = type
c no DF if dftype = 'none' or 'Jean'
        dist( icmp ) = ( idftyp( icmp ) .gt. 1 ) .and.
     +                 ( idftyp( icmp ) .ne. 28 )
        if( dist( icmp ) .and. ( .not. term ) .and. master )write( no,
     +               '( '' DF type selected is '', a4, '' index'', i3 )'
     +                                           )type, idftyp( icmp )
c get DF parameters
        j = npars( idftyp( icmp ) )
        dfcns( 1, icmp ) = 0
        dfcns( 2, icmp ) = 0
        if( .not. term )then
c need initial Q at least for a disc
          if( disc( icmp ) .and. ( .not. dist( icmp ) ) )j = 1
          if( j .gt. 0 )then
            do while ( line( m:m ) .ne. ' ' )
              m = m + 1
              if( m .gt. n )go to 1
            end do
            read( line( m:n ), *, err = 1 )
     +                                    ( dfcns( i, icmp ), i = 1, j )
            if( disc( icmp ) .and. idftyp( icmp ) .eq. 1 )then
              if( master )write( no,
     +              '( '' Initial Q of disc'', f7.4 )' )dfcns( 1, icmp )
            else
              if( master )write( no,
     +    '( '' DF constants'', 2f9.4 )' )( dfcns( i, icmp ), i = 1, j )
            end if
          end if
        end if
        retract = .false.
        if( type .eq. 'ZANG' )then
c input Q for a Mestel disc
          if( term )call gtreal( 'Enter nominal Q', dfcns( 1, icmp ) )
          if( cmpmas( icmp ) .gt. 0. )then
            dfcns( 1, icmp ) = 7. /
     +                      ( cmpmas( icmp ) * dfcns( 1, icmp ) )**2 - 1
          else
            dfcns( 1, icmp ) = 7. / dfcns( 1, icmp )**2 - 1
          end if
        else if( j .gt. 0 )then
          if( term )then
            if( j .eq. 1 )call gtreal(
     +                           'Enter DF constant', dfcns( 1, icmp ) )
            if( j .eq. 2 )call gtreals(
     +                   'Enter two DF constants', dfcns( 1, icmp ), 2 )
          end if
        else
          if( disc( icmp ) .and. ( .not. dist( icmp ) ) )then
            if( term )call gtreal( 'Enter initial Q', dfcns( 1, icmp ) )
            if( dfcns( 1, icmp ) .lt. 0. )call crash( 'MSETUP',
     +                                      'Nonsense Q value entered' )
            kcold = dfcns( 1, icmp ) .eq. 0
          end if
        end if
c set tapers for some disks
        if( disc( icmp ) )then
          Lztapr( icmp ) = .false.
          if( cdft( icmp ) .ne. 'CPB0' )Lzrmn( icmp ) = 0
c power law inner angular momentum taper for some discs
          if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +        ( cdft( icmp ) .eq. 'EVAN' ) )then
            if( term )then
              call gtreal(
     +        'Enter index for inner angular momentum taper, suggest 4',
     +                                                  indexi( icmp ) )
            else
              call getline( ni, line )
              if( line( 1:4 ) .ne. 'IIND' )then
                if( master )print *, 'Inner taper index line missing?'
                go to 1
              end if
              n = lnblnk( line )
              read( line( 11:n ), *, err = 1 )indexi( icmp )
            end if
            if( master )print *, 'Inner taper index is',indexi( icmp )
            Lztapr( icmp ) = .true.
          else
c Lz_crit for flipping orbits in other discs except for Miyamoto & Omega models
            if( .not. ( ( cdft( icmp ) .eq. 'MIYA' ) .or.
     +                  ( cdft( icmp ) .eq. 'OMEG' ) .or.
     +                  ( cdft( icmp ) .eq. 'CPB0' ) ) )then
c Kalnajs DF for the isochrone disk is a special case
              if( ( ctype( icmp ) .eq. 'ISOC' ) .and.
     +            ( cdft( icmp ) .eq. 'KALN' ) )then
                if( term )then
                  retract = gtlogl(
     +                      'Do you want Kalnajs rule for retro stars' )
                else
                  call getline( ni, line )
                  if( line( 1:4 ) .ne. 'KRET' )then
                    if( master )print *,
     +               'Kret line for Kalnajs retro rule logical missing?'
                    go to 1
                  end if
                  n = lnblnk( line )
                  read( line( 11:n ), *, err = 1 )retract
                end if
              end if
c other disks or DFs
              if( Lzcrt( icmp ) .eq. 0. )then
                Lzrmn( icmp ) = 0
                lskip = ( ( cdft( icmp ) .eq. 'NONE' ) .or.
     +                    ( cdft( icmp ) .eq. 'EVCO' ) ) .or.
     +                  ( ( ctype( 1 ) .eq. 'POWE' ) .and.
     +                    ( dfcns( 3, icmp ) .ge. 0. ) ) .or. retract
                if( .not. lskip )then
                  if( term )then
                    call gtdble( 'Enter Lz_crit', Lzcrt( icmp ) )
                  else
                    call getline( ni, line )
                    if( line( 1:4 ) .ne. 'LZCR' )then
                      if( master )print *, 'Lzcrit line missing'
                      go to 1
                    end if
                    n = lnblnk( line )
                    read( line( 11:n ), *, err = 1 )Lzcrt( icmp )
                  end if
c turn on taper if this is activated
                  Lztapr( icmp ) = Lzcrt( icmp ) .gt. 0.
                end if
              end if
            end if
          end if
c outer angular momentum taper for truncated models
          if( trnctd( icmp ) )then
            if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +          ( cdft( icmp ) .eq. 'EVAN' ) .or.
     +          ( cdft( icmp ) .eq. 'EVCO' ) )then
c power law tapers for power law models
              if( term )then
                call gtreal(
     +       'Enter index for outer angular momentum taper (suggest 6)',
     +                                                  indexo( icmp ) )
                call gtdble(
     +              'Enter mean radius for outer angular momentum taper'
     +                                  // ' (suggest rmax - 5)', rtap )
              else
                call getline( ni, line )
                if( line( 1:4 ) .ne. 'OUTT' )go to 1
                n = lnblnk( line )
                read( line( 11:n ), *, err = 1 )indexo( icmp ), rtap
              end if
              Lztapr( icmp ) = .true.
c
              if( cdft( icmp ) .eq. 'EVCO' )then
                if( gtlogl(
     +       'Do you want to remove high energy radial orbits?' ) )then
                  indexi( icmp ) = 1
                else
                  indexi( icmp ) = 0
                end if
              end if
            else
c simple cubic taper
              if( term )then
                call gtdble(
     +               'Enter radius of inner edge of outer taper', rtap )
                Lztapr( icmp ) = sngl( rtap ) .lt. rtrunc( icmp )
              else
                call getline( ni, line )
                if( line( 1:4 ) .ne. 'TAPE' )go to 1
                n = lnblnk( line )
                read( line( 11:n ), *, err = 1 )Lztapr( icmp )
                if( Lztapr( icmp ) )read(
     +                    line( 11:n ), *, err = 1 )Lztapr( icmp ), rtap
              end if
            end if
          end if
        end if
c active mass fraction
        fmass( icmp ) = 1
        if( trnctd( icmp ) )then
          if( ctype( icmp ) .eq. 'UNKN' )then
            fmass( icmp ) = 1
          else
            r2 = rtrunc( icmp )
            if( cmpmas( icmp ) .gt. 0. )then
              if( disc( icmp ) )then
c too soon to integrate over DF so force analytic expression
                local = dist( icmp )
                dist( icmp ) = .false.
                fmass( icmp ) = gmassd( r2 ) / cmpmas( icmp )
                dist( icmp ) = local
              else
                if( ctype( icmp ) .eq. 'ASDI' )then
                  j = ncmp + 1 - icmp
                  fmass( icmp ) = fmass( j )
                else
                  if( ctype( icmp ) .eq. 'ADIA' )then
                    fmass( icmp ) = 1
                  else
                    fmass( icmp ) = gmassh( r2 ) / cmpmas( icmp )
                  end if
                end if
              end if
            end if
          end if
        end if
      end do
c initialize compressed model
      if( jcomp .gt. 0 )then
        icmp = jcomp
        call compin( .true. )
        cc = cc0( icomp )
        dfcns( 3, icmp ) = cc
        cmppar( 1, icmp ) = cc
      end if
c check for nonsense parameters for an ASDI halo
      if( froz )then
        j = 0
        m = 0
        do i = 1, ncmp
          if( ctype( i ) .ne. 'ASDI' )then
            if( disc( i ) )then
              if( ( m .gt. 0 ) .and.
     +             ( ctype( i ) .ne. ctype( m ) ) )call crash( 'MSETUP',
     +                  'ASDI halo for more than 1 distinct component' )
              m = i
            else
              call crash( 'MSETUP',
     +      'ASDI halo for a model with some other non-disc component' )
            end if
          end if
        end do
      end if
c set variables in / model /
      call inimod
c set outer boundary
      cc = 0
      do i = 1, ncmp
        cc = max( cc, rtrunc( i ) )
      end do
      call cutoff( cc )
      icmp = 0
c set disk taper values - need vcirc for whole model
      do icmp = 1, ncmp
        if( disc( icmp ) )then
          if( cdft( icmp ) .eq. 'SHUE' )dist( icmp ) = .false.
          if( Lzcrt( icmp ) .gt. 0. )then
            Lzrmn( icmp ) = -min( rmax * vcirc( rmax ),
     +                                      Lzcrt( icmp ) )
          end if
c outer angular momentum taper for truncated models
          if( trnctd( icmp ) )then
            if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +          ( cdft( icmp ) .eq. 'EVAN' ) .or.
     +          ( cdft( icmp ) .eq. 'EVCO' ) )then
              Lzmno( icmp ) = rtap * vcirc( rtap )
              if( master )print *, 'Outer taper index is',
     +                     indexo( icmp ), ' centered on', Lzmno( icmp )
            else
c simple cubic taper
              r2 = rtrunc( icmp )
              if( dist( icmp ) )then
                Lztmx( icmp ) = r2 * vcirc( r2 )
                Lztmn( icmp ) = Lztmx( icmp )
                if( Lztapr( icmp ) )then
                  Lztmn( icmp ) = rtap * vcirc( rtap )
                  if( master )print
     +                 '( '' Cubic taper in angular momentum from'', '//
     +           'f10.4, '' to'', f10.4 )', Lztmn( icmp ), Lztmx( icmp )
                end if
              else
                Lztmx( icmp ) = r2
                Lztmn( icmp ) = Lztmx( icmp )
                if( Lztapr( icmp ) )then
                  Lztmn( icmp ) = rtap
                  if( master )print
     +                 '( '' Cubic taper in radius from'', '//
     +           'f10.4, '' to'', f10.4 )', Lztmn( icmp ), Lztmx( icmp )
                end if
              end if
            end if
          end if
          if( retract )Lzrmn( icmp ) = -Lztmx( icmp )
        end if
c revise fmass
        if( ( cmpmas( icmp ) .gt. 0. ) .and. Lztapr( icmp ) .and.
     +      ( ctype( icmp ) .ne. 'UNKN' ) )then
          fmass( icmp ) = gmassi( r2 ) / cmpmas( icmp )
        end if
        if( cdft( icmp ) .eq. 'SHUE' )dist( icmp ) = .true.
      end do
      icmp = 0
      return
c error reading data card
    1 if( master )then
        print '( 2x, a )', line( 1:60 )
        print *,
     + 'This error could simply mean an expected input line was missing'
      end if
      call crash( 'MSETUP', 'Problem with input data line' )
      stop
c intentional abort
    2 call crash( 'MSETUP', 'END string input' )
      stop
      end
