      subroutine dfhead( unit, irec, out, warm )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to write or read the header of a .dfn file.
c The warm parameter is ignored for writes, but for reads:
c   if( warm )then
c     the DF is assumed to be for mass component icmp in common / model /
c     and most of the header values are ignored
c     this is the normal call from psetup
c   else
c     the values read are placed in common blocks and some are reported out
c     this call is used from dflook etc
c   end if
c
c calling arguments
      integer irec, unit
      logical out, warm
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
      character*4 dftype, dtype, htype
      logical gtlogl
      real version
      real*8 splint2, vcirc
c
c local arrays
      integer jdf( 4 ), jdum( 3 )
      logical done( mcmp )
      real dum( 17 )
      save done
      equivalence ( dum( 1 ), jdum( 1 ) )
c
c local variables
      integer ia, ib, ip, j, k, norb, np
      logical flag, la, ldsc, lnum, luqm
      real a, indi, indo, ret, rmx, v
      real*8 Lzmn, Lzmo, Lzmx, r1, r2
      equivalence ( la, ia ), ( ldsc, ib )
c
      data done / mcmp * .false. /
c
      if( no .eq. 0 )no = 6
      v = version( 0. )
c
      rewind unit
      if( out )then
        if( ( icmp .lt. 1 ) .or. ( icmp .gt. ncmp ) )then
          print *, icmp
          call crash( 'DFHEAD', 'Unrecognized component number' )
        end if
        write( no, * )'Writing distribution function header to unit',
     +                 unit
        ret = Lzcrt( icmp )
        if( retract )ret = -1
        rmx = rmax
        write( unit )irec, ( jdfcs( j, icmp ), j = 1, 4 ), ncmp,
     +               rmx, ret, rhole, Lztmn( icmp ), Lztmx( icmp ),
     +               Lzmno( icmp ), indexi( icmp ), indexo( icmp ),
     +               icmp, ncode, uqmass, cvers, numcal, dum
        do ip = 1, ncmp
          write( unit )imtyp( ip ), impar( ip ), idftyp( ip ),
     +                 ( dfcns( j, ip ), j = 1, 3 ), cmpmas( ip ),
     +                 fmass( ip ), rscale( ip ),
     +                 rtrunc( ip ) / rscale( ip ), disc( ip )
        end do
      else
        if( master )write( no, * )
     +            'Reading distribution function header from unit', unit
c read header records
        if( warm )then
c read minimal information
          read( unit, err = 1 )irec, ( jdf( j ), j = 1, 4 ), np,
     +             rmx, ret, rhole, Lzmn, Lzmx, Lzmo, indi, indo, ip, j,
     +             luqm, v, numcal
          ip = 0
          do while ( .not. done( icmp ) )
            ip = ip + 1
            if( ip .gt. np )call crash( 'DFHEAD',
     +                      'No matching component found in .dfn file' )
            read( unit )( dum( j ), j = 1, 10 ), ldsc
c allow for differing binary representations of .true.!
            la = .true.
            if( ldsc )ib = ia
c identify population that matches info expected
            flag = .true.
            flag = flag .and. ( jdum( 1 ) .eq. imtyp( icmp ) )
            flag = flag .and. ( jdum( 2 ) .eq. impar( icmp ) )
            flag = flag .and. ( jdum( 3 ) .eq. idftyp( icmp ) )
c skip if these parameters have been encountered before for another component
            if( flag .and. ( icmp .gt. 1 ) )then
              if( idftyp( icmp ) .eq. idftyp( ip ) )flag =
     +                                   flag .and. ( .not. done( ip ) )
            end if
c check other properties
            if( flag )then
              a = dum( 4 ) - dfcns( 1, icmp )
              if( dum( 4 ) .ne. 0. )a = a / dum( 4 )
              flag = flag .and. abs( a ) .lt. 1.e-5
c              flag = flag .and. ( dum( 5 ) .eq. dfcns( 2, icmp ) )
c              flag = flag .and. ( dum( 6 ) .eq. dfcns( 3, icmp ) )
c              flag = flag .and. ( dum( 7 ) .eq. cmpmas( icmp ) )
c              flag = flag .and. ( dum( 9 ) .eq. rscale( icmp ) )
              dum( 10 ) = dum( 10 ) * rscale( icmp )
              flag = flag .and.
     +                  ( abs( dum( 10 ) - rtrunc( icmp ) ) .lt. 1.e-5 )
              flag = flag .and. ( ldsc .eqv. disc( icmp ) )
              done( icmp ) = .true.
            end if
            fmass( icmp ) = dum( 8 )
            done( icmp ) = done( icmp ) .or. ( ip .eq. np )
          end do
          if( .not. flag )then
            if( master )then
              print *, 'DF info inconsistent with pop', icmp
              print *, ( jdum( j ), j = 1, 3 )
              print *, imtyp( icmp ), impar( icmp ), idftyp( icmp )
              print *, ( dum( j ), j = 4, 6 )
              print *, ( dfcns( j, icmp ), j = 1, 3 )
              print *, dum( 7 ), dum( 9 ), dum( 10 ), ldsc
              print *, cmpmas( icmp ), rscale( icmp ), rtrunc( icmp ),
     +                 disc( icmp )
            end if
            flag = gtlogl( 'Is this difference acceptable' )
            if( .not. flag )call crash( 'DFHEAD', 'Wrong .dfn file' )
          end if
          rmx = max( rmx, dum( 10 ) )
c skip any unread header records
          if( ip .lt. np )then
            do j = ip + 1, np
              read( unit )
            end do
          end if
        else
c read and interpret full header information
          read( unit, err = 1 )irec, ( jdf( j ), j = 1, 4 ), ncmp,
     +               rmx, ret, rhole, Lzmn, Lzmx, Lzmo, indi, indo,
     +               icmp, ncode, luqm, v, numcal, dum
          if( ( icmp .lt. 1 ) .or. ( icmp .gt. mcmp ) )then
            if( master )print *, 'icmp =', icmp
            call crash( 'DFHEAD', 'Invalid component number' )
          end if
          k = -1
          do ip = 1, ncmp
            read( unit, err = 1 )imtyp( ip ), impar( ip ), idftyp( ip ),
     +                       ( dfcns( j, ip ), j = 1, 3 ), cmpmas( ip ),
     +                          fmass( ip ), rscale( ip ), rtrunc( ip ),
     +                          disc( ip )
            if( imtyp( ip ) .eq. 23 )k = ip
c initialize character strings
            j = imtyp( ip )
            if( disc( ip ) )then
              ctype( ip ) = dtype( j )
              if( ( ctype( ip ) .eq. 'SOFT' ) .or.
     +            ( ctype( ip ) .eq. 'GAUS' ) )epsiln = dfcns( 3, ip )
            else
              ctype( ip ) = htype( j )
            end if
            cdft( ip ) = dftype( idftyp( ip ) )
            rtrunc( ip ) = rscale( ip ) * rtrunc( ip )
          end do
c initialize uncompressed halo and start over
          if( k .gt. 0 )then
            ip = icmp
            icmp = k
            call compin( .true. )
            icmp = ip
            rewind unit
            read( unit )irec, ( jdf( j ), j = 1, 4 ), ncmp,
     +                  rmx, ret, rhole, Lzmn, Lzmx, Lzmo, indi, indo,
     +                  icmp, ncode, luqm, v, dum
            do ip = 1, ncmp
              read( unit )imtyp( ip ), impar( ip ), idftyp( ip ),
     +                    ( dfcns( j, ip ), j = 1, 3 ), cmpmas( ip ),
     +                    fmass( ip ), rscale( ip ), rtrunc( ip ),
     +                    disc( ip )
            end do
c select component
            do ip = 1, ncomp
              k = iccmp( ip )
              if( k .eq. icmp )icomp = ip
            end do
          end if
          lnum = numcal
        end if
        if( uqmass .neqv. luqm )then
          if( master )print *, 'warning uqmass reset to', luqm
          uqmass = luqm
        end if
c assign values for selected component
        norb = 1
        do j = 1, 4
          jdfcs( j, icmp ) = jdf( j )
          if( j .lt. 4 )norb = norb * jdf( j )
        end do
        Lztmn( icmp ) = Lzmn
        Lztmx( icmp ) = Lzmx
        Lzmno( icmp ) = Lzmo
        indexi( icmp ) = indi
        indexo( icmp ) = indo
c interpret ret
        retract = ret .lt. 0.
        if( retract )then
          Lzcrt( icmp ) = 0
        else
          Lzcrt( icmp ) = ret
        end if
        if( warm )then
c check outer boundary only
          if( ( .not. cmprssd ) .and.
     +        ( abs( rmx - rtrunc( icmp ) ) .gt. .1 ) )then
            if( master )print *, 'Warning truncation radius of comp',
     +  icmp, ' is', rtrunc( icmp ), ' but rmax read in DFHEAD is', rmx
            if(
     +        gtlogl( 'Accept the value in the .dfn file header' ) )then
              rtrunc( icmp ) = rmx
              rmx = 0
              do ip = 1, ncmp
                rmx = max( rmx, rtrunc( ip ) )
              end do
              call cutoff( rmx )
            end if
          end if
        else
c initialize variables in / model /
          call inimod
c set truncation and angular momentum tapers
          Lztapr( icmp ) = ( Lzcrt( icmp ) .gt. 0.d0 ) .or.
     +                     ( ( Lztmn( icmp ) .gt. 0.d0 ) .and.
     +                       ( Lztmn( icmp ) .lt. Lztmx( icmp ) ) )
        end if
c get table of numerical potential, if it exists and is not already read
        if( ( cdft( icmp ) .eq. 'SHUE' ) .and. numcal .and.
     +      ( .not. lgfld ) )then
          k = -1
          call opnfil( k, 'ftb', 'unformatted', 'old', 'seq', ip )
          if( ip .eq. 0 )then
            if( .not. allocated( arad ) )allocate ( arad( mradad ) )
            if( .not. allocated( Phin ) )allocate ( Phin( mradad ) )
            if( .not. allocated( frin ) )allocate ( frin( mradad ) )
            read( k )nradad, arad, Phin, frin
            close( k )
c initialize splines
            if( .not. allocated( plamda ) )then
              allocate ( plamda( nradad + 4 ) )
              allocate ( phinc( nradad + 4 ) )
            end if
            if( .not. allocated( flamda ) )then
              allocate ( flamda( nradad + 4 ) )
              allocate ( frinc( nradad + 4 ) )
            end if
            r2 = arad( nradad / 2 )
            r1 =
     +          splint2( arad, Phin, nradad, r2, plamda, phinc, .true. )
            r1 =
     +          splint2( arad, frin, nradad, r2, flamda, frinc, .true. )
            lgfld = .true.
          else
            call crash( 'DFHEAD', 'numerical potential file not found' )
          end if
        end if
        call cutoff( rmx )
        if( retract )then
          Lzrmn( icmp ) = -Lztmx( icmp )
        else
          Lzrmn( icmp ) = 0
          if( Lzcrt( icmp ) .gt. 0.d0 )Lzrmn( icmp ) =
     +                       -min( rmax * vcirc( rmax ), Lzcrt( icmp ) )
        end if
        if( master )then
          write( no, * )norb, ' particles'
          if( jdfcs( 3, icmp ) .gt. 1 )write( no, * )jdfcs( 3, icmp ),
     +                                      ' particles for each E & Lz'
          write( no,
     +             '( 5x, '' Mass represented by particles ='', f12.6 )'
     +                                                    )fmass( icmp )
c inner angular momentum taper
          if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +        ( cdft( icmp ) .eq. 'EVAN' ) )then
c power law tapers for power law models
            write( no, * )' Inner taper index is ', indexi( icmp )
          else
            write( no, '( '' h_crit ='', f8.3 )' )Lzcrt( icmp )
          end if
c outer angular momentum taper
          if( ( cdft( icmp ) .eq. 'ZANG' ) .or.
     +        ( cdft( icmp ) .eq. 'EVAN' ) .or.
     +        ( cdft( icmp ) .eq. 'EVCO' ) )then
            write( no, * )' Outer taper index is ', indexo( icmp ),
     +                                    ' centred on ', Lzmno( icmp )
          else
c simple cubic taper
            if( Lztmn( icmp ) .lt. Lztmx( icmp ) )then
              write( no, '( '' Cubic taper in angular momentum from'','
     +       // ' f10.4, '' to'', f10.4 )' )Lztmn( icmp ), Lztmx( icmp )
            else
              write( no, '( '' Sharp cut-off in angular momentum at'','
     +                                      // ' f10.4 )' )Lztmx( icmp )
            end if
          end if
        end if
      end if
      return
c header records could not be read
    1 print *, 'Could this be an old .dfn file?'
      call crash( 'DFHEAD',  'Failed to read header' )
      stop
      end
