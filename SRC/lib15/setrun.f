      subroutine setrun
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to read the parameters from the .dat file needed to set up a run
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
c externals
      character*4 htype
      integer lnblnk
      logical gtlogl
      real roundup
      real*8 vcirc
c
c local arrays
      integer inst, minst, ninst
      parameter ( ninst = 26, minst = 7, inst = minst * mcmp )
      character*4 a( ninst + 1 )
      logical ldist( mcmp ), set( ninst ), setp( minst, mcmp )
      logical lhev( mgrids )
c
c local variables
      character line*80
      integer i, icard, i1, i2, j, k, norbit !, nactiv
      logical ld, lh
      real actmas, b, gm, Lzmax, u
      real*8 r
      include 'inc/pi.f'
c
      data a / 'LSCA', 'ICMP', 'TIME', 'ANAL', 'SAVE', 'UQMA', 'ZONE',
     +         'GUAR', 'CNTD', 'PERT', 'SUPP', 'PLOT', 'LGSP', 'DVEL',
     +         'HVEL', 'HBIN', 'SPHB', 'ZANL', 'MONI', 'ZPRF', 'RNGS',
     +         'RHOR', 'OFFG', 'REST', 'FRQS', 'END ', 'trap' /
      data set / ninst * .false. /, setp / inst * .false. /
c
      if( ncode .gt. 11 )call crash( 'SETRUN',
     +                                        'Unrecognized grid type' )
c initialise / reslts /
      do i = 1, ndattyp
        dattyp( i ) = .false.
      end do
      lsave = .false.
c RVLF time centering
      lfrv = .true.
c set default nprop
      if( twod )nprop = 5
      if( threed )nprop = 9
c nsect should be set already
      if( nsect .le. 0 )call crash( 'SETRUN', 'nsect not set' )
c default is not to recenter
      netcom = .false.
      netlmo = .false.
c default is not to randomize particles during an unload
      lurand = .false.
c
c read and interpret next data card
      if( master)write( no,
     +                '( / '' Reading model parameters in SETRUN'' )'  )
      line( 1:4 ) = '    '
      do while ( line( 1:3 ) .ne. 'END' )
        call getline( ni, line )
        icard = 1
        do while ( line( 1:4 ) .ne. a( icard ) )
          icard = icard + 1
c unrecognised code word
          if( icard .gt. ninst )then
            if( master )then
              print '( / '' Keyword is '', a4 )', line( 1:4 )
              write( no, '( / ''Keyword is '', a4 )' )line( 1:4 )
            end if
            call crash( 'SETRUN', 'Unrecognised data card' )
            line( 1:4 ) = a( ninst + 1 )
          end if
        end do
        if( icard .le. ninst )then
          if( ( icard .ne. 2 ) .and. set( icard ) .and. master )print
     + '( '' Warning: repeat data card '', a4, '//
     +                       ' '' encountered in SETRUN'' )', a( icard )
          set( icard ) = .true.
        end if
c
c scale length in grid units
        if( line( 1:4 ) .eq. 'LSCA' )then
          read( line( 11:80 ), *, err = 1 )lscale
c information about population
        else if( line( 1:4 ) .eq. 'ICMP' )then
          read( line( 11:80 ), *, err = 1 )icmp
          if( ( icmp .lt. 1 ) .or. ( icmp .gt. mcmp )
     +                        )call crash( 'SETRUN', 'Impossible ICMP' )
          do while ( line( 1:3 ) .ne. 'END' )
            call getline( ni, line )
c number of particles in this population
            if( line( 1:4 ) .eq. 'NPAR' )then
              read( line( 11:80 ), *, err = 1 )nsp( icmp )
              rigidp( icmp ) = nsp( icmp ) .eq. 0
              testp( icmp ) = cmpmas( icmp ) .eq. 0.
              npr( icmp ) = 0
              setp( 1, icmp ) = .true.
c start type instructions for given population
            else if( line( 1:4 ) .eq. 'STAR' )then
              ldist( icmp ) = dist( icmp )
              read( line( 11:80 ), *, err = 1 )
     +                          dist( icmp ), smr( icmp ), quiet( icmp )
              setp( 2, icmp ) = .true.
c grid number for this population
            else if( line( 1:4 ) .eq. 'PGRD' )then
              read( line( 11:80 ), *, err = 1 )igrd( icmp )
              setp( 3, icmp ) = .true.
c vertical thickness of population
            else if( line( 1:4 ) .eq. 'Z0IN' )then
              read( line( 11:80 ), *, iostat = i )z0init( icmp ),
     +                                            iztyp( icmp )
              if( i .ne. 0 )then
                read( line( 11:80 ), *, err = 1 )z0init( icmp )
c Gaussian vertical distribution is default - ignored for non-discs
                iztyp( icmp ) = 2
              end if
              setp( 4, icmp ) = .true.
c population center
            else if( line( 1:4 ) .eq. 'CPOP' )then
              read( line( 11:80 ), *, err = 1 )
     +                                     ( comi( i, icmp ), i = 1, 3 )
              setp( 5, icmp ) = .true.
c population velocity
            else if( line( 1:4 ) .eq. 'VPOP' )then
              read( line( 11:80 ), *, err = 1 )
     +                                     ( comi( i, icmp ), i = 4, 6 )
              setp( 6, icmp ) = .true.
c optional disk outer taper
            else if( line( 1:4 ) .eq. 'TAPE' )then
c taper only for selected component
              Lztapr( icmp ) = .false.
              if( cdft( icmp ) .ne. 'CPB0' )Lzrmn( icmp ) = 0
c outer angular momentum taper for truncated models
              if( disc( icmp ) .and. trnctd( icmp ) )then
c simple cubic taper
                r = rtrunc( icmp )
                Lztmx( icmp ) = r * vcirc( r )
                Lztmn( icmp ) = Lztmx( icmp )
                read( line( 11:80 ), *, err = 1 )r
                if( sngl( r ) .lt. rtrunc( icmp ) )then
                  Lztmn( icmp ) = r * vcirc( r )
                  Lztapr( icmp ) = .true.
                  if( master )print
     +                 '( '' Cubic taper in angular momentum from'', '//
     +           'f10.4, '' to'', f10.4 )', Lztmn( icmp ), Lztmx( icmp )
                else
                  if( master )print *, 'inner taper from', sngl( r ),
     +                       ' but disk is truncated at', rtrunc( icmp )
                  call crash( 'SETRUN', 'Nonsense outer taper radius' )
                end if
              end if
              if( .not. Lztapr( icmp ) )then
                if( master )write( no, * )'Population', icmp
                call crash( 'SETRUN',
     +                              'Outer taper not sensible for pop' )
              end if
              setp( 7, icmp ) = .true.
            end if
          end do
          line( 1:4 ) = 'ICMP'
c length of timestep
        else if( line( 1:4 ) .eq. 'TIME' )then
          read( line( 11:80 ), *, err = 1 )ts
c iphys determines the number of time-steps between model analyses
        else if( line( 1:4 ) .eq. 'ANAL' )then
          read( line( 11:80 ), *, err = 1 )iphys
c physical data are to be saved on a file
        else if( line( 1:4 ) .eq. 'SAVE' )then
          lsave = .true.
c read instructions for results file
          do while ( line( 1:4 ) .ne. 'END ' )
            call getline( ni, line )
            if( line( 1:4 ) .eq. 'DANL' )danl = .true.
            if( line( 1:4 ) .eq. 'DNST' )dnst = .true.
            if( line( 1:4 ) .eq. 'FRQS'
     +             )call crash( 'SETRUN', 'FRQS now needs a parameter' )
c            if( line( 1:4 ) .eq. 'GREY' )grey = .true.
            if( line( 1:4 ) .eq. 'INTG' )intg = .true.
            if( line( 1:4 ) .eq. 'MOIT' )moit = .true.
            if( line( 1:4 ) .eq. 'PNTL' )pntl = .true.
            if( line( 1:4 ) .eq. 'S3DC' )s3dc = .true.
            if( line( 1:4 ) .eq. 'LVAL' )lval = .true.
          end do
          line( 1:4 ) = a( 6 )
c variable particle masses
        else if( line( 1:4 ) .eq. 'UQMA' )then
          read( line( 11:80 ), *, err = 1 )uqmass
c number and radii of zones
        else if( line( 1:4 ) .eq. 'ZONE' )then
          read( line( 11:80 ), *, err = 1 )nzones
          nstep( 1 ) = 1
          tzone( 1 ) = 1.
          if( nzones .gt. 1 )then
c read time-step factor and radius limit for each zone
            do i = 2, nzones
              call getline( ni, line )
              read( line( 11:80 ), *, err = 1 )nstep( i ), zbr( i )
              tzone( i ) = real( nstep( i )**2 )
            end do
          end if
c guard radii
        else if( line( 1:4 ) .eq. 'GUAR' )then
          read( line( 11:80 ), *, err = 1 )nguard
          if( ( nguard .lt. 0 ) .or. ( nguard .gt. mguard ) )then
            if( master )then
              write( no, * )'Impossible value for nguard', nguard
              print *, 'Impossible value for nguard', nguard
            end if
            call crash( 'SETRUN', 'Guard radii' )
          end if
          if( nguard .eq. 0 )then
            rguard = 0
          else
c read time-step factor and radius limit for guard region
            do i = 1, nguard
              call getline( ni, line )
              read( line( 11:80 ), *, err = 1 )nsub( i ), gbr( i )
              if( nsub( i ) .eq. 0 )then
                if( master )then
                  print '( 3x, a )', line( 11:80 )
                  print *,
     +                'data on guard radii in old order - please revise'
                end if
                read( line( 11:80 ), *, err = 1 )gbr( i ), nsub( i )
              end if
            end do
            rguard = gbr( 1 )
          end if
c grid recentering option
        else if( line( 1:4 ) .eq. 'CNTD' )then
          read( line( 11:80 ), *, iostat = i )icstp, j
c centroid rule is the default
          if( i .ne. 0 )then
            read( line( 11:80 ), *, err = 1 )icstp
            j = 1
          end if
c set flag for centering and time interval
          centrd = .true.
          icstp = max( icstp, 1 )
          if( ( ngrid .eq. 1 ) .and.
     +        ( noslfg .or. dr3d .or. bht ) )call crash( 'SETRUN',
     +                    'Grid centering requested when none exists!' )
c determine which rule
          ltroid = j .eq. 1
          lbind = j .eq. 2
          lheavy = j .eq. 3
          if( ltroid )then
c k = 1/2 option only for now
            kci = 1
          else if( lbind )then
c assume 100 particles to define center
            nebind = 100
          else if( .not. lheavy )then
            print *, icstp, j
            call crash( 'SETRUN', 'Unknown centering rule' )
          end if
c perturbing or central mass option
        else if( line( 1:4 ) .eq. 'PERT' )then
          pertbn = .true.
          i = 11
          do while ( line( i:i ) .eq. ' ' )
            i = i + 1
          end do
          extbar = line( i:i+2 ) .eq. 'BAR'
          exspir = line( i:i+3 ) .eq. 'SPIR'
          exgnrc = line( i:i+3 ) .eq. 'GNRC'
          if( .not. ( extbar .or. exspir .or. exgnrc ) )then
            isphp = -1
            do j = 1, nhalos
              if( line( i:i+3 ) .eq. htype( j ) )isphp = j
            end do
            exsphr = isphp .gt. 0
            if( .not. exsphr )call crash( 'SETRUN',
     +                                        'Unrecognized perturber' )
          end if
          if( exsphr .or. exgnrc )then
            i = i + 4
            if( exgnrc )then
              read( line( i:80 ), *, err = 1 )mptbr
            else
              read( line( i:80 ), *, err = 1 )mptbr, eptbr
            end if
c initial position and velocity
            call getline( ni, line )
            if( line( 1:4 ) .ne. 'POSI' )call crash( 'SETRUN',
     +                       'Initial position of perturber not found' )
            read( line( 11:80 ), *, err = 1 )xptrb
            call getline( ni, line )
            if( line( 1:4 ) .ne. 'VELO' )call crash( 'SETRUN',
     +                       'Initial velocity of perturber not found' )
            read( line( 11:80 ), *, err = 1 )vptrb
          else if( extbar )then
            i = i + 3
            if( quapp )then
              read( line( i:80 ), *, err = 1 )mptbr, abar, omegab,
     +                                                        bbar, cbar
              eptbr = bbar / abar
c set alpha2 & beta2
              call qudini
            else
              read( line( i:80 ), *, err = 1 )mptbr, abar, omegab, eptbr
              bbar = abar * eptbr
            end if
            bphase = 0
          else if( exspir )then
            i = i + 4
            read( line( i:80 ), *, err = 1 )omegsp, epspi
            sphase = 0
          else
            call crash( 'SETRUN', 'Unrecognized perturber' )
          end if
c correction (supplementary) forces for whole model
        else if( line( 1:4 ) .eq. 'SUPP' )then
          read( line( 11:80 ), *, err = 1 )suppl
          lsupst = .false.
c recentering options at outset only
        else if( line( 1:4 ) .eq. 'REST' )then
          read( line( 11:80 ), *, err = 1 )netcom, netlmo
c step interval for plot file
        else if( line( 1:4 ) .eq. 'PLOT' )then
          read( line( 11:80 ), *, err = 1 )ipict
          plot = ipict .gt. 0
c radial spacing - number of values per lscale
        else if( line( 1:4 ) .eq. 'FRQS' )then
          read( line( 11:80 ), *, err = 1 )i
          frqs = .true.
          drfrqs = lscale / real( i )
c logarithmic spirals to be determined - check space first
        else if( line( 1:4 ) .eq. 'LGSP' )then
          read( line( 11:80 ), *, err = 1 )i1, i2
c ensure i1 is odd
          i1 = 2 * ( ( i1 - 1 ) / 2 ) + 1
          np = i1
          nm = i2
          lgsp = np .gt. 0
c control parameters for analysis of disc particle velocity field
        else if( line( 1:4 ) .eq. 'DVEL' )then
          read( line( 11:80 ), *, err = 1 )i1, i2
          vfld = i1 + i2 .gt. 0
          ndrad = max( i1, 1 )
          ndang = max( i2, 1 )
c control parameters for analysis of halo particle velocity field
        else if( line( 1:4 ) .eq. 'HVEL' )then
          read( line( 11:80 ), *, err = 1 )i1, i2
          vflh = i1 + i2 .gt. 0
          nhrbin = max( i1, 1 )
          nhzbin = max( i2, 1 )
c ensure an odd number of vertical bins
          nhzbin = 2 * ( ( nhzbin - 1 ) / 2 ) + 1
c number of angular momentum bins
        else if( line( 1:4 ) .eq. 'HBIN' )then
          read( line( 11:80 ), *, err = 1 )i1
          angm = i1 .gt. 0
          maxLz = max( i1 - 2, 1 )
c spherical Bessel functions
        else if( line( 1:4 ) .eq. 'SPHB' )then
          read( line( 11:80 ), *, err = 1 )i1, i2
          jnmax = i1
          jlmax = i2
          if( s3d .or. scf )then
            jlmax = min( jlmax, s3lmax )
          else
            if( jlmax .gt. s3maxl )then
              if( master )then
                print *, jlmax, ' requested', s3maxl, ' declared'
                write( no, * )jlmax, ' requested', s3maxl, ' declared'
              end if
              call crash( 'SETRUN', 's3maxl in / grids / too small' )
            end if
            s3lmax = jlmax
            s3ntm = ( jlmax + 1 ) * ( jlmax + 1 ) / 2
          end if
          sphb = jnmax .gt. 0
c azimuthal mid-plane variations - check space first
        else if( line( 1:4 ) .eq. 'ZANL' )then
          read( line( 11:80 ), *, err = 1 )i1, i2
          nrz = i1
          nmz = max( i2, 1 )
          zanl = nrz .gt. 0
c particle monitoring for integral conservation - input is stride
        else if( line( 1:4 ) .eq. 'MONI' )then
          read( line( 11:80 ), *, err = 1 )i1
          moni = i1 .gt. 0
          nmskip = i1
c vertical density profile analysis for disk
        else if( line( 1:4 ) .eq. 'ZPRF' )then
          read( line( 11:80 ), *, err = 1 )i1, i2
          zprf = i1 + i2 .gt. 0
c ensure an odd number of vertical bins
          i2 = 2 * ( ( i2 - 1 ) / 2 ) + 1
c check space
          i1 = max( i1, 1 )
          nbzr = max( i1, 1 )
          nbzz = i2
c rings option
        else if( line( 1:4 ) .eq. 'RNGS' )then
          rngs = .true.
          read( line( 11:80 ), *, iostat = k )nrings, i, j
          if( k .ne. 0 )then
            if( ngrid .eq. 1 )then
              read( line( 11:80 ), *, err = 1 )nrings, i
              j = 1
            else
              call crash( 'SETRUN', 'Must specify grid for rings' )
            end if
          end if
          ncmp = ncmp + 1
          if( ncmp .gt. mcmp )call crash( 'SETRUN', 'Invalid pop #' )
          npring = i
          nsp( ncmp ) = nrings * npring
          igrd( ncmp ) = j
          z0init( ncmp ) = 0
          iztyp( ncmp ) = 2
          do i = 1, minst
            setp( i, ncmp ) = .true.
          end do
c density profile option
        else if( line( 1:4 ) .eq. 'RHOR' )then
          rhor = .true.
          read( line( 11:80 ), *, err = 1 )mprho
          mprho = max( mprho, 1 )
c off-grid particle options
        else if( line( 1:4 ) .eq. 'OFFG' )then
          read( line( 11:80 ), *, err = 1 )offanal, offfor, offthf,
     +                                     offmnp
c          read( line( 11:80 ), *, iostat = i )offanal, offfor, offthf,
c     +                                        offmnp
c          if( i .ne. 1 )then
c            offmnp = .false.
c            read( line( 11:80 ), *, iostat = i )offanal, offfor, offthf
c            if( i .ne. 1 )then
c              offthf = .false.
c              read( line( 11:80 ), *, iostat = i )offanal, offfor
c              if( i .ne. 1 )then
c                offfor = .false.
c                read( line( 11:80 ), *, err = 1 )offanal
c              end if
c            end if
c          end if
        end if
      end do
c
c data stream need not set values for every option
      if( master )write( no, * )
      ld = .false.
      lh = .false.
!      nactiv = 0
      do icmp = 1, ncmp
c set look up table for active components
!        if( nsp( icmp ) .gt. 0 )then
!          nactiv = nactiv + 1
!          jactiv( icmp ) = nactiv
!        end do
        if( disc( icmp ) )then
          ld = .true.
        else
          lh = .true.
        end if
c check consistency of set up
        if( quiet( icmp ) )then
          if( netcom .or. netlmo )call crash( 'SETRUN',
     +                 'Arrest model inconsistent with a quiet start' )
          set( 24 ) = .true.
        end if
        if( rigidp( icmp ) .and. ( netcom .or. netlmo ) )write( no,
     +        * )'Arrest model requested with a rigid component present'
      end do
c report options that were not set
      do i = 6, ninst - 1
        if( .not. set( i ) )then
          if( i .eq. 6 )then
            if( master )write( no,
     +                  * )'Particles assumed to all have the same mass'
            uqmass = .false.
          else if( i .eq. 7 )then
            if( master )write( no, * )'No time step zones set'
            nzones = 1
            nstep( 1 ) = 1
            tzone( 1 ) = 1
          else if( i .eq. 8 )then
            if( master )write( no, * )'No guard radii set'
            nguard = 0
            rguard = 0
          else if( i .eq. 9 )then
            if( master )write( no, * )'Grid center fixed'
            centrd = .false.
            icstp = 1
          else if( i .eq. 10 )then
            if( master )write( no, * )'No external perturber'
            pertbn = .false.
          else if( i .eq. 11 )then
            if( master )write( no, * )'No supplementary forces'
            suppl = .false.
            lsupst = .false.
          else if( i .eq. 12 )then
            if( master )write( no,
     +                           * )'Particle data for pictures omitted'
            plot = .false.
          else if( i .eq. 13 )then
            if( ld .and. master )write( no,
     +                          * )'Logarithmic spiral analysis omitted'
            lgsp = .false.
          else if( i .eq. 14 )then
            if( ld .and. master )write( no,
     +                         * )'Disk velocity field analysis omitted'
            vfld = .false.
          else if( i .eq. 15 )then
            if( lh .and. master )write( no,
     +                         * )'Halo velocity field analysis omitted'
            vflh = .false.
          else if( i .eq. 16 )then
            if( master )write( no, * )'Angular momentum binning omitted'
            angm = .false.
          else if( i .eq. 17 )then
            if( lh .and. master )write( no,
     +                   * )'Spherical Bessel function analysis omitted'
            sphb = .false.
          else if( i .eq. 18 )then
            if( ld .and. master )write( no, *
     +                                  )'Disk bending analysis omitted'
            zanl = .false.
          else if( i .eq. 19 )then
            if( master )write( no,
     +                           * )'Particle monitoring option omitted'
            moni = .false.
          else if( i .eq. 20 )then
            if( ld .and. master )write( no, *
     +                                  )'Disk z-profile option omitted'
            zprf = .false.
          else if( i .eq. 21 )then
            if( master )write( no, * )'No rings'
            rngs = .false.
          else if( i .eq. 22 )then
            if( master )write( no, * )'No density profile'
            rhor = .false.
          else if( i .eq. 23 )then
            if( master )write( no, * )'Off grid particles discarded'
            offanal = .false.
            offfor = .false.
            offthf = .false.
            offmnp = .false.
          else if( i .eq. 24 )then
            if( master )write( no,
     +                         * )'Initial model will not be recentered'
          end if
          set( i ) = .true.
        end if
      end do
      do icmp = 1, ncmp
c if only 1 grid, the population must be on it!
        if( ngrid .eq. 1 )then
          igrd( icmp ) = 1
          setp( 3, icmp ) = .true.
        end if
c z thickness meaningless for 2-D or spherical components
        if( ( ( .not. disc( icmp ) ) .and.
     +        ( .not. sphrod( icmp ) ) ) .or. twod )then
          setp( 4, icmp ) = .true.
          iztyp( icmp ) = 0
          if( twod )then
            z0init( icmp ) = 0
          else
            z0init( icmp ) = 1
          end if
        end if
c start type, grid number, initial thickness and inital Q need not be set
c   for a rigid mass distribution
        if( rigidp( icmp ) )then
          setp( 2, icmp ) = .true.
          setp( 3, icmp ) = .true.
          setp( 4, icmp ) = .true.
        end if
        do i = 1, minst
          if( .not. setp( i, icmp ) )then
            if( master )then
              if( i .eq. 1 )then
                write( no, * )'Number of particles not set'
                print *, 'Number of particles not set'
              else if( i .eq. 2 )then
                write( no, * )'Start type not set'
                print *, 'Start type not set'
              else if( i .eq. 3 )then
                write( no, * )'Grid number not set'
                print *, 'Grid number not set'
              else if( i .eq. 4 )then
                write( no, * )'Initial thickness not set'
                print *, 'Initial thickness not set'
              end if
            end if
            if( i .eq. 5 )then
              if( master )write( no, * )'Center of mass at grid centre'
              do j = 1, 3
                comi( j, icmp ) = 0.
              end do
            else if( i .eq. 6 )then
              if( master )write( no, * )'No initial net momentum'
              do j = 4, 6
                comi( j, icmp ) = 0.
              end do
            else if( i .eq. 7 )then
              ld = disc( icmp ) .and. trnctd( icmp )
              if( ld .and. master )then
                if( Lztapr( icmp ) )then
                  write( no, * )'Outer disc density tapered'
                else
                  write( no, * )'Sharp disc truncation'
                end if
              end if
            end if
            if( i .lt. 5 )then
              if( master )then
                write( no, '( ''Population'', i3 )' )icmp
                print *, 'Population', icmp
              end if
              call crash( 'SETRUN', 'Insufficent data for population' )
            end if
          else
            if( i .eq. 5 )then
              if( master )write( no,
     +        '( 9x, ''Centre for population'', i3, '' is'', 3f10.3 )' )
     +                               icmp, ( comi( j, icmp ), j = 1, 3 )
            else if( i .eq. 6 )then
              if( master )write( no,
     +   '( ''Initial velocity of population'', i3, '' is'', 3f10.3 )' )
     +                               icmp, ( comi( j, icmp ), j = 4, 6 )
            end if
          end if
        end do
      end do
c check whether all options set
      ld = .false.
      do i = 1, ninst
        if( .not. set( i ) )ld = .true.
      end do
      if( ld )then
        if( master )then
          do i = 1, ninst
            if( .not. set( i ) )write( no,
     +                                '( 3x, a4, '' missing'' )' )a( i )
          end do
        end if
        call crash( 'SETRUN', 'Insufficient data' )
      end if
c radial extent of rings
      if( rngs )then
        rtrunc( ncmp ) = 0
        do i = 1, ncmp - 1
          rtrunc( ncmp ) = max( rtrunc( i ), rtrunc( ncmp ) )
        end do
        cmpmas( ncmp ) = 0
        npr( ncmp ) = npring
        imtyp( ncmp ) = 19
        disc( ncmp ) = .true.
        z0init( ncmp ) = 0
        dfcns( 1, ncmp ) = 0
        dist( ncmp ) = .false.
        smr( ncmp ) = .false.
        quiet( ncmp ) = .true.
        rscale( ncmp ) = 1
        testp( ncmp ) = .true.
        rigidp( ncmp ) = .false.
      end if
c report rigid mass components
      if( master )then
        write( no, * )
        do icmp = 1, ncmp
          if( rigidp( icmp ) )then
            if( disc( icmp ) )then
              write( no, * )'rigid disc selected for pop', icmp
            else
              write( no, * )'rigid halo selected for pop', icmp
            end if
          end if
        end do
c print out perturber info if appropriate
        if( pertbn )then
          if( exsphr )then
            if( exsphr )write( no,
     +                       * )'Spherical perturber ' // htype( isphp )
            write( no, '( ''Perturber mass and core radius'', 2f8.3 )' )
     +               mptbr, eptbr
            write( no, '( '' Initial position of perturber'', 3f8.3 )' )
     +               xptrb
          else if( extbar )then
            write( no,
     +    '( ''Bar mass, semi-major axis and pattern speed'', 3f8.3 )' )
     +               mptbr, abar, omegab
          else if( exspir )then
            write( no,
     +  '( ''Spiral pertbn with pattern speed and strength'', 2f8.3 )' )
     +              omegsp, epspi
          else
            if( .not. exgnrc
     +                 )call crash( 'SETRUN', 'Unrecognized perturber' )
          end if
        end if
      end if
c set offset pointer for heavy particles
      call switch( 0 )
      if( dr3d )then
c ensure drctN particles are at the top if other methods are in use too
        kdrct = 0
        if( ngrid .gt. 1 )then
          i = igrd( ncmp )
          if( igrid( i ) .ne. 10 )call
     +          crash( 'SETRUN', 'Heavy particles not last population' )
          do i = 1, ncmp - 1
            kdrct = kdrct + nsp( i )
          end do
        end if
      end if
c total number of particles
      nbod = 0
      actmas = 0
      do i = 1, ncmp
        nbod = nbod + nsp( i )
        if( nsp( i ) .gt. 0 )actmas = actmas + cmpmas( i )
      end do
      if( actmas .le. 0. )then
        print *, 'No active mass'
        if( .not. gtlogl( 'Do you want to continue' ) )call
     +                               crash( 'SETRUN', 'No active mass' )
        actmas = 1
      end if
c logic for heavy particles if needed
      if( centrd .and. lheavy )then
        if( .not. uqmass )call crash( 'SETRUN',
     +  'Heavy particle option for grid centering requires uqmass = T' )
        do i = 1, ngrid
c ensure only one per grid
          lhev( i ) = .true.
          do j = 1, ncmp
            if( lhev( i ) .and. ( igrd( j ) .eq. i ) )then
              lhev( i ) = .false.
c remember which component was selected
              jebind( i ) = j
              nsp( j ) = nsp( j ) + 1
              nbod = nbod + 1
            end if
          end do
        end do
      end if
c set default grid size for gridless codes
      if( sf2d .or. sf3d .or. scf )then
        minmaxr = minmaxr * lscale * rtrunc( 1 )
        maxmaxr = maxmaxr * lscale * rtrunc( 1 )
        rgrid( 1 ) = max( rgrid( 1 ), maxmaxr )
      end if
c convert and check zoning information
      rguard = rguard * lscale
      zbr( 1 ) = max( rguard, rinh( 1 ) )
      if( nzones .gt. 1 )then
c RVLF required if zones are activated
        lfrv = .true.
        do i = 2, nzones
          zbr( i ) = zbr( i ) * lscale
          if( zbr( i ) .lt. zbr( i - 1 ) )then
            call crash( 'SETRUN', 'Zone boundaries in illogical order' )
          end if
          if( zbr( i ) .lt. rinh( 1 ) )then
            if( master )then
              print *, 'Zone', i
              write( no, * )'Zone', i
            end if
            call crash( 'SETRUN', 'boundary inside guard radius' )
          end if
          if( zbr( i ) .gt. rgrid( ngrid ) )then
            if( master )then
              print *, 'Zone', i, zbr( i ), rgrid( ngrid )
              write( no, * )'Zone', i, zbr( i ), rgrid( ngrid )
            end if
            call crash( 'SETRUN', 'boundary outside grid' )
          end if
        end do
      end if
      mzone = 1
c convert and check guard radii information
      if( nguard .gt. 0 )then
        do i = 1, nguard
          gbr( i ) = gbr( i ) * lscale
          if( i .gt. 1 )then
            if( gbr( i ) .gt. gbr( i - 1 ) )
     +     call crash( 'SETRUN', 'Guard boundaries in illogical order' )
          end if
        end do
      else
        if( rinh( 1 ) .gt. 0. )call crash( 'SETRUN',
     +                    'Must set step length for particles in hole' )
      end if
      if( ( nzones .gt. 1 ) .or. ( nguard .gt. 0 ) )then
        sphbnd = noslfg .or. s3d .or. scf .or. p3d
        if( ( .not. sphbnd ) .and.
     +        .not. ( twod .or. p3a .or. sf3d ) )call crash( 'SETRUN',
     +                                   'Unknown zone boundary shape' )
      end if
c compute sizes of arrays ans set some defaults
      call setptr
c check that mass fraction is consistent with particle number for equal mass
      icmp = ncmp
      i = nbod
      do k = 1, ncmp
        if( testp( k ) )then
          icmp = icmp - 1
          i = i - nsp( k )
        end if
      end do
      if( ( icmp .gt. 1 ) .and. ( .not. uqmass ) )then
        do j = 1, icmp
          if( .not. rigidp( j ) )then
            k = real( i ) * cmpmas( j ) / actmas + .5
            if( nsp( j ) .ne. k )then
              if( master )then
                print *,
     +                 'Inconsistent number of particles in population',
     +                                           j, ' for mass fraction'
                print *, 'nsp =', nsp( j ), ' required =', k,
     +                          ' mass fraction =', cmpmas( j ) / actmas
              end if
            end if
          end if
        end do
      end if
c prevent attempts to analyse the model at substeps of the maximum
      if( mod( iphys, nstep( nzones ) ) .ne. 0 )call crash( 'SETRUN',
     +           'Analysis not at a multiple of largest particle step' )
c point plotting for pictures
      if( plot )ipict = max( ipict, iphys )
      zshift = 0
      do i = 1, ncmp
        u = 2. * rtrunc( i ) * lscale
        zshift = max( zshift, u )
      end do
      zshift = max( zshift, 2. * zm( ngrid ) )
c obsolete orbit tracking option
      norbit = 0

c scale model such that:
c         G is 1, length scale is 1 and time unit is 1
      gm = actmas * lscale**3 * ts**2
c pmass = scaled mass of 1 particle - ignore massless ring particles
      i = nbod
      if( rngs )i = nbod - nsp( ncmp )
      if( i .eq. 0 )then
        if( noslfg )then
c set notional value
          pmass = 1
        else
          call crash( 'SETRUN', 'No massive particles' )
        end if
      else
        pmass = gm / real( i )
      end if
c gvfac = conversion factor for velocities
      gvfac = 1. / ( lscale * ts )
c output type of leap-frog selected
      if( master )then
        if( lfrv )then
          write( no, * )'Retarded velocity leap-frog scheme selected'
        else
          write( no, * )'Advanced position leap-frog scheme selected'
        end if
      end if
c zoning
      if( master )then
        if( nzones .gt. 1 )then
          write( no,
     + '( '' Zoning instructions:  beyond radius  step factor'', ' //
     +  '10( / 22x, f10.3, i10 ) )' )
     +    ( zbr( i ), nstep( i ), i = 2, nzones )
        else
          write( no, * )'All particles have the same time step'
        end if
        if( nguard .gt. 0 )then
          write( no,
     +    '( ''Steps subdivided when r<'', f10.4 )' )rguard / lscale
          do i = 1, nguard
            write( no, '( i12, '' substeps inside'', f10.4 )'
     +                                  )nsub( i ), gbr( i ) / lscale
          end do
        end if
      end if
c centroiding
      if( centrd )then
        icstp = max( icstp, nstep( nzones ) )
        if( master )write( no,
     +       '( ''Centroid of particle distribution revised every'' ' //
     +                                       ', i4, '' steps'' )' )icstp
      end if
c grid interpolation scheme
      if( master )then
        if( c3d )then
          if( jmass .eq. 1 )write( line( 1:4 ), '( a )' )' NGP'
          if( jmass .eq. 2 )write( line( 1:4 ), '( a )' )' CIC'
          if( jmass .eq. 3 )write( line( 1:4 ), '( a )' )' TSC'
          write( line( 5:27 ), '( a )' )' interpolation selected'
          write( no, '( a )' )line( 1:27 )
        else if( p2d .or. p3d .or. p3a .or. s3d .or. sf3d )then
          if( jmass
     +             .eq. 1 )write( line( 1:13 ), '( a )' )'  First order'
          if( jmass
     +             .eq. 2 )write( line( 1:13 ), '( a )' )' Second order'
          if( jmass
     +             .eq. 3 )write( line( 1:13 ), '( a )' )'     Accurate'
          write( line( 14:36 ), '( a )' )' interpolation selected'
          write( no, '( a )' )line( 1:36 )
        else
          write( no, * )'No relevant interpolation scheme for this code'
        end if
      end if
c off-grid particle options
      if( offanal )then
        ld = .false.
        if( offfor )then
          if( master )write( no, * )
     +             'Off-grid particles accelerated by the direct method'
          ld = offthf .or. offmnp
        else if( offthf )then
          if( master )write( no, * )
     +     'Off-grid particles accelerated using theoretical expression'
          ld = offmnp
        else if( offmnp )then
          if( master )write( no, * )
     +   'Off-grid particles accelerated using point mass approximation'
        else
          if( master )write( no, * )
     +           'Off-grid particles advanced without being accelerated'
        end if
        if( ld )call crash( 'SETRUN',
     + 'Multiple acceleration methods selected for off-grid particles' )
      end if
c scaled units
      do icmp = 1, ncmp
        if( master )write( no,
     +               '( '' Population'', i3, '' mass is'', f5.2, ' //
     +               ' '' represented by'', i9, '' particles'' )' )icmp,
     +               cmpmas( icmp ), nsp( icmp )
        if( threed .and. disc( icmp ) )then
          if( master )write( no,
     +     '( '' Initial scale height of disc'', f6.2 )' )z0init( icmp )
          if( ( z0init( icmp ) .gt. 0. ) .and.
     +        ( nsp( icmp ) .gt. 0 )  .and.
     +  ( ( iztyp( icmp ) .lt. 1 ) .or. ( iztyp( icmp ) .gt. 5 ) ) )then
            if( master )then
              print *, 'iztyp =', iztyp( icmp )
              write( no, * ) 'iztyp =', iztyp( icmp )
            end if
            call crash( 'SETRUN', 'Unrecognized iztype' )
          end if
          if( master )then
            if( iztyp( icmp ) .eq. 1 )write( no, * )
     +                              ' Spitzer analytic isothermal sheet'
            if( iztyp( icmp ) .eq. 2 )write( no, * )
     +                                      ' Gaussian vertical profile'
            if( iztyp( icmp ) .eq. 3 )write( no, * )
     +            ' sech^2 vertical profile, but not Spitzer dispersion'
            if( iztyp( icmp ) .eq. 4 )write( no, * )
     +                           ' rounded exponential vertical profile'
            if( iztyp( icmp ) .eq. 5 )write( no, * )
     +                                   ' exponential vertical profile'
          end if
          znorm( icmp ) = 1
        end if
c call out a possible inconsistency in .dat file
        if( ldist( icmp ) .neqv. dist( icmp ) )then
          if( master )then
            print *, 'mismatch for DF of component', icmp
            print *, 'model information has DF', ldist( icmp ),
     +              ' while setrun read DF', dist( icmp )
          end if
          if( disc( icmp ) )then
            if( .not. gtlogl( 'Do you want to continue' ) )call crash(
     +                                      'SETRUN', 'Abort selected' )
c reset cutoff parameters for disk, if no DF is preferred
            if( .not. dist( icmp ) )call tapset
          end if
        end if
        if( dist( icmp ) .and. master )write( no,
     +'( '' Particle coordinates expected on file run'', i4, ''.dfn'' )'
     +                                                             )irun
      end do
      if( master )then
        if( suppl )write( no, * )
     +         'Radial forces from particles corrected to analytic form'
        if( nsect .gt. 1 )write( no,
     +                       '( i6, ''-fold symmetry imposed'' )' )nsect
        write( no,
     +         '( '' Length scale is'', f7.3, '' grid units'' )' )lscale
        write( no, '( 11x, '' Time step is'', f8.4 )' )ts
      end if
c compute self-energy of a possible rigid halo
      call halset
      if( ( norbit .ne. 0 ) .and. master )write( no,
     +         '( / i8, '' particle orbits will be followed'' )' )norbit
      norbit = 2 * norbit
c prepare for spherical Bessel analysis
      if( sphb )then
        call sphbjz( jlmax )
        rbess = rgrid( ngrid )
      end if
c particle velocity field analysis
      if( vfld )then
        i = ndrad
      else if( vflh )then
        i = nhrbin
      else
c best to set a notional value
        i = 1
      end if
      drfac = real( i ) / lscale
      if( vfld )then
c outer edge for disc analysis
        u = -1
        do i = 1, ncmp
          if( disc( i ) .and. ( nsp( i ) .gt. 0 ) )then
            b = 2 * lscale * rtrunc( i )
            if( b .gt. u )then
              u = b
              j = igrd( i )
            end if
          end if
        end do
        if( u .le. 0.
     +    )call crash( 'SETRUN', 'Disc analysis for no live component' )
        u = min( u, rgrid( j ) )
        ndrad = u * drfac
      end if
      if( vflh )then
c outer edge for halo analysis
        u = -1
        do i = 1, ncmp
          if( ( .not. disc( i ) ) .and. ( nsp( i ) .gt. 0 ) )then
            b = 1.2 * lscale * rtrunc( i )
            if( b .gt. u )then
              u = b
              j = igrd( i )
            end if
          end if
        end do
        if( u .le. 0.
     +    )call crash( 'SETRUN', 'Halo analysis for no live component' )
        u = min( u, rgrid( j ) )
        nhrbin = u * drfac
      end if
      if( vfld )then
        thfac = real( ndang ) / ( 2. * pi )
        idskp = ndang * nprop
      end if
      if( vflh )then
        i = real( nhzbin ) / real( nhrbin )
        hzfac = .5 * real( i + 1 )
        hzfac = drfac * hzfac
        ihskp = nhzbin * nprop
      end if
c frequency analysis - set a default scale if none was set above
      if( .not. frqs )drfrqs = lscale
c angular momentum bins
      if( angm )then
        Lzmax = 0
        do i = 1, ncmp
          if( nsp( i ) .gt. 0 )Lzmax = max( Lzmax, 1.2 * rtrunc( i ) )
        end do
        Lzmax = roundup( Lzmax )
        Lzmax = Lzmax * lscale**2 * ts
        Lzfac = real( maxLz ) / Lzmax
      end if
c z-density profile to go 2 radial bins beyond initial disk edge
      if( zprf )then
        i = max( 1, nbzr -2 )
        rbzf = real( i ) / ( rtrunc( 1 ) * lscale )
        zbzf = real( nbzz / 2 ) / zm( 1 )
      end if
c particle monitoring
      if( moni )then
        i = ( nbod - 1 ) / nmskip + 1
        nmonit = i
        nmskip = nmskip * nwpp
      end if
      if( lsave .and. master )then
        write( no, '( '' Analysis every'', i4, '' steps to save:'' )' )
     +                  iphys
        if( dnst )write( no, * )'      density array'
        if( danl )write( no, * )'      density transform'
        if( s3dc )write( no, * )'      coeffs of s3d grid'
        if( grey )write( no, * )'      grey scale array'
        if( pntl )write( no, * )'      potential array'
        if( frqs )write( no,
     +'( 7x, ''frequencies at radial spacing'', f8.3 )' )drfrqs / lscale
        if( moit )write( no, * )'      moment of inertia tensors'
        if( orbs )write( no, * )'      selected orbits'
        if( lgsp )write( no,
     + '( 7x, ''Logarithmic spiral analysis into'', i4, ' //
     + '  '' Fourier harmonics and'', i4, '' pitch angles'' )' )nm, np
        if( sphb )write( no,
     + '( 7x, ''Spherical Bessel analysis up to order'', i4, ' //
     + '  '' with'', i4, ' //
     + '  '' radial functions, and radius scale'', f8.3 )' )jlmax,
     +         jnmax, rbess / lscale
        if( intg )write( no, * )'      global integrals'
        if( moni )write( no, '( 7x, ''integrals for'', i8,'//
     +                       ' '' particles'' )' )nmonit
        if( angm )write( no, * )'      angular momentum distribution'
        if( vfld )write( no,
     + '( 7x, ''Disc velocity field analysis into'', i4, ' //
     + '  '' radial and'', i4, '' azimuthal bins'', 16x, ' //
     + ' ''with outer boundary at'', f6.2, '' scales'' )' )ndrad, ndang,
     +                                real( ndrad ) / ( lscale * drfac )
        if( vflh )write( no,
     + '( 7x, ''Halo velocity field analysis into'', i4, ' //
     + ' '' radial and'',  i4, '' vertical bins'', 10x / ' //
     + ' ''with radial and vertical boundaries at'', 2f7.3, ' //
     + ' '' scales'' )' )nhrbin, nhzbin,
     +     real( nhrbin + 1 ) / ( lscale * drfac ),
     +    ( real( nhzbin / 2 ) + .5 ) / ( lscale * hzfac )
        if( rhor )write( no, '( 7x, ''halo density from every'', i6,' //
     +                       ' ''th particle'' )'  )mprho
        if( zprf )write( no,
     + '( 7x, ''z-density profile analysis into'', i4, ' //
     + '  '' radial and'', i4, '' vertical bins'' )' )nbzr, nbzz
        if( plot )then
          write( no, '( '' Points for pictures every'', i4, ' //
     +                                            ' '' steps'' )' )ipict
        end if
      end if
c set up self-energy and self-force terms
      call slfset
c initialise / stfile /
      noff = 0
      do i = 1, 4
        ncen( i ) = 0
      end do
      nshort = 0
c initialise off angular momentum
      do i = 1, 2
        do j = 1, 3
          angoff( j, i ) = 0
          amoff( j, i ) = 0
        end do
      end do
c sensible defaults
      istep = 0
      phys = .true.
c flag successful completion
      if( master )then
        write( no, * )'SETRUN successfully completed'
        write( no, * )
      end if
      return
    1 if( master )print *, 'Problem with data card:'
      i = lnblnk( line )
      if( master )print *, line( 1:i )
      call crash( 'SETRUN', 'Error reading data' )
      stop
      end
