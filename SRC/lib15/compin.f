      subroutine compin( in )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to read the information in the pre-calculated compress.dat file
c   for an adiabatically compressed halo, to scale the model and set up the
c   spline fits
c
c calling argument
      logical in
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      character*4 dftype, dtype, htype
      integer nchars
      logical gtlogl
      real*8 distfn, E0tab, Phihal, splint2
c
c local arrays
      character*4 sdft( 2 ), stype( 2 )
      integer itemp( 2, 2 ), jmtyp( 2 )
      logical ldsc( 2 )
      real temp( 7, mcmp )
c
c local variables
      character fname*20
      integer i, idisc0, ih, iuse, j, jcmp, k, kcmp, n, unit
      logical flag
      real hmass, hscale, ratio, sl, sm, sr
      real*8 extmas, extscl, r, x
      save iuse
c
      data iuse / -1 /
c
      if( in )then
        irun = max( irun, 0 )
        if( iuse .ne. irun )then
c select file name to use
          if( ( icmp .lt. 1 ) .or. ( icmp .gt. ncmp )
     +                )call crash( 'COMPIN', 'Nonsense value for icmp' )
c open data file
          if( irun .gt. 0 )then
            call opnfil( unit, 'cmp', 'unformatted', 'old', 'seq', i )
            if( i .ne. 0 )call crash( 'COMPIN', '.cmp file not found' )
          else
            ih = imtyp( icmp )
            if( ih .eq. 27 )then
              fname = 'VKA1.dat'
            else
              fname = 'compress.dat'
            end if
            n = nchars( fname )
            open( unit, file = fname( 1:n ), form = 'unformatted',
     +            status = 'old', iostat = i )
            if( i .ne. 0 )call crash( 'COMPIN',
     +                                   fname(1:n)//' file not found' )
          end if
c read header records
          read( unit, iostat = j )nradad, ncomp, nrigid, ihalo0, idfn0,
     +                     iccmp, ircmp, cc0, dfm0, dfb0, r0core, r1cusp
          if( j .eq. 0 )then
            kcmp = ncomp + nrigid
            read( unit )( ( temp( j, i ), j = 1, 4 ), i = 1, kcmp ),
     +                          ( ldsc( i ), jmtyp( i ), i = 1, nrigid )
            if( master )then
              print *
              print *, 'Header information read by compin:'
              print *, 'nradad, ncomp, nrigid', nradad, ncomp, nrigid
              print *, 'r0core, r1cusp', r0core, r1cusp
              do jcmp = 1, kcmp
                print *, 'Component', jcmp
                do i = 1, ncomp
                  if( iccmp( i ) .eq. jcmp )then
                    j = ihalo0( i )
                    k = idfn0( i )
                    print '( 1x, a, 2( i3, 1x, a ) )',
     +                     'is a compressed halo, original type and DF',
     +                                     j, htype( j ), k, dftype( k )
                    print '( 1x, a, 3f10.4 )',
     +  'and parameters cc0, dfm0, dfb0', cc0( i ), dfm0( i ), dfb0( i )
                  end if
                end do
                do i = 1, nrigid
                  if( ircmp( i ) .eq. jcmp )then
                    j = jmtyp( i )
                    if( ldsc( i ) )then
                      print *, 'is a rigid disk', j, dtype( j )
                    else
                      print *, 'is a rigid sphere', j, htype( j )
                    end if
                  end if
                end do
                print '( 1x, a, 4f10.4 )',
     +       'rtrunc, mass, scale, fmass', ( temp( j, jcmp ), j = 1, 4 )
              end do
              print *
            end if
          else

            if( master )print *, 'reading old header'
            call crash( 'COMPIN', 'read old header is untested code' )

c attempt to read older header
            rewind unit
            read( unit, iostat = j )nradad, idisc0, ihalo0( 1 ),
     +                  idfn0( 1 ), cc0( 1 ), dfm0( 1 ), dfb0( 1 ),
     +                  r0core, r1cusp, extscl, extmas,
     +                  temp( 1, icmp ), hmass, hscale
            if( j .ne. 0 )then
c very old header
              rewind unit
              read( unit, iostat = j )nradad, idisc0, ihalo0( 1 ),
     +                  idfn0( 1 ), cc0( 1 ), dfm0( 1 ), dfb0( 1 ),
     +                  r0core, r1cusp, extscl, extmas,
     +                  temp( 1, icmp )
              if(
     +         j .ne. 0 )call crash( 'COMPIN', 'Failed to read header' )
c old defaults
              hscale = 1
              hmass = 1
            end if
c supply missing data
            ncomp = 1
            nrigid = 1
            if( ncmp .eq. 2 )then
c choose which is the disk
              j = 1
              if( icmp .eq. 1 )j = 2
            else
              if( master )then
                print *, 'model has', ncmp, ' components'
                print *, 'compressible halo is component', icmp
              end if
              call crash( 'COMPIN', 'I do not know what to do' )
            end if
c set rigid disk
            ircmp( 1 ) = j
            jmtyp( j ) = idisc0
            ldsc( j ) = .true.
            temp( 1, j ) = rtrunc( j )
            temp( 2, j ) = extmas
            temp( 3, j ) = extscl
c set compressed halo
            iccmp( 1 ) = icmp
            temp( 1, icmp ) = rtrunc( icmp )
            temp( 2, icmp ) = hmass
            temp( 3, icmp ) = hscale
          end if
c check values and parameters from the compress.dat file
          if( kcmp .ne. ncmp )then
            if( master )then
              print *, 'Your current model has', ncmp, ' components'
              print *, 'while the compress.dat file has', kcmp
            end if
            call crash( 'COMPIN', 'Inconsistent number of cmps' )
          end if
          flag = .true.
          do i = 1, nrigid
            j = ircmp( i )
            flag = flag .and. ( imtyp( j ) .eq. jmtyp( i ) )
            flag = flag .and. ( disc( j ) .eqv. ldsc( i ) )
          end do
          do i = 1, ncomp
            j = iccmp( i )
            flag = flag .and. ( imtyp( j ) .eq. 23 )
            flag = flag .and. ( idftyp( j ) .eq. 25 )
          end do
          if( flag )then
          else
c try to determine whether re-ordering is all that is needed
            if( ncmp .eq. 2 )then
              flag = .true.
              j = 3 - iccmp( 1 )
              flag = flag .and. ( imtyp( j ) .eq. 23 )
              flag = flag .and. ( idftyp( j ) .eq. 25 )
              j = 3 - ircmp( 1 )
              flag = flag .and. ( imtyp( j ) .eq. jmtyp( 1 ) )
              flag = flag .and. ( disc( j ) .eqv. ldsc( 1 ) )
              if( flag )then
                if( master )print *, 'swap is enough'
c swap components in compress data to be consistent with input .dat file
                iccmp( 1 ) = 3 - iccmp( 1 )
                ircmp( 1 ) = 3 - ircmp( 1 )
                do k = 1, 4
                  x = temp( k, 1 )
                  temp( k, 1 ) = temp( k, 2 )
                  temp( k, 2 ) = x
                end do
              else
                if( master )then
                  j = 3 - iccmp( 1 )
                  print *, j, imtyp( j ), 23
                  print *, j, idftyp( j ), 25
                  j = 3 - ircmp( 1 )
                  print *, j, imtyp( j ), jmtyp( 1 )
                  print *, j,  disc( j ), ldsc( 1 )
                  print *,
     +              'compress.dat model does not match your input model'
                end if
                call crash( 'COMPIN', 'Fix it' )
              end if
            else
              if( master )then
                print *, 'The', ncmp, ' components in this model differ'
     +  // ' between the compress.dat file'
                print *, ' and other input information, perhaps just'
     +  // ' in their order'
              end if
              call crash( 'COMPIN', 'fix it' )
            end if
          end if
c read in compressed halo arrays
          if( nradad .gt. mradad )call
     +                      crash( 'COMPIN', 'Common arrays too small' )
          allocate ( arad( mradad ) )
          allocate ( phin( mradad ) )
          allocate ( gmn( mradad, ncomp ) )
          allocate ( rhon( mradad, ncomp ) )
          read( unit )( arad( i ), i = 1, nradad ),
     +                ( phin( i ), i = 1, nradad ),
     +                ( ( gmn( i, j ), i = 1, nradad ), j = 1, ncomp ),
     +                ( ( rhon( i, j ), i = 1, nradad ), j = 1, ncomp )
          close( unit )
c set rigid component(s) - not always set
          do k = 1, nrigid
            j = ircmp( k )
            rigidp( j ) = .true.
            imtyp( j ) = jmtyp( k )
            disc( j ) = ldsc( k )
            if( disc( j ) )then
              ctype( j ) = dtype( imtyp( j ) )
            else
              ctype( j ) = htype( imtyp( j ) )
            end if
c use the value read from file if fmass is not preset
            if( fmass( j ) .le. 0. )fmass( j ) = temp( 4, j )
          end do
c not everything is read from headers for compressible components
          do k = 1, ncomp
            j = iccmp( k )
            rigidp( j ) = .false.
            disc( j ) = .false.
            dist( j ) = .true.
            ctype( j ) = 'ADIA'
            cdft( j ) = 'COMP'
c use the value read from file if fmass is not preset
            if( fmass( j ) .le. 0. )fmass( j ) = temp( 4, j )
          end do
c accept values from input .dat file, but check for consistent scaling
          flag = .true.
          ratio = cmpmas( 1 ) / temp( 2, 1 )
          do i = 2, ncmp
            flag = flag .and.
     +                abs( ratio * temp( 2, i ) - cmpmas( i ) ) .lt. .01
          end do
          ratio = rscale( 1 ) / temp( 3, 1 )
          do i = 2, ncmp
            flag = flag .and.
     +                abs( ratio * temp( 3, i ) - rscale( i ) ) .lt. .01
          end do
          do i = 1, ncmp
c            flag = flag .and.
c     +                abs( ratio * temp( 1, i ) - rtrunc( i ) ) .lt. .01
            if( abs( ratio * temp( 1, i ) - rtrunc( i ) ) .gt. .01 )then
              if( master )then
                print *, 'Warning: truncation radius discrepancy for' //
     +                  ' mass component', i
                print *, temp( 1, i ), rtrunc( i )
              end if
              rtrunc( i ) = temp( 1, i )
            end if
          end do
          if( .not. flag )then
            if( master )then
              print *,
     +    'compressed model cannot be scaled to parameters in .dat file'
              do i = 1, ncmp
                print *, 'parameters for component', i
                print *, 'truncation radii', rtrunc( i ), temp( 1, i )
                print *, '   radial scales', rscale( i ), temp( 3, i )
                print *, '          masses', cmpmas( i ), temp( 2, i )
              end do
            end if
            call crash( 'COMPIN', 'inconsistent model' )
          end if
c rescale arrays
          j = -1
          do i = 1, ncmp
            if( disc( i ) )then
              if( j .gt. 0 )call crash( 'COMPIN',
     +                                      'Multiple disk components' )
              j = i
            end if
          end do
          if( j .gt. 0 )then
            sl = rscale( j ) / temp( 3, j )
            sm = cmpmas( j ) / temp( 2, j )
            sr = sm / sl**3
            if( ( abs( sl - 1. ) .gt. .001 ) .or.
     +          ( abs( sm - 1. ) .gt. .001 ) .or.
     +          ( abs( sr - 1. ) .gt. .001 ) )then
              do i = 1, nradad
                arad( i ) = arad( i ) * sl
                phin( i ) = phin( i ) * sm / sl
                do icomp = 1, ncomp
                  gmn( i, icomp ) = gmn( i, icomp ) * sm
                  if( r1cusp )then
                    rhon( i, icomp ) = rhon( i, icomp ) + log10( sr )
                  else
                    rhon( i, icomp ) = rhon( i, icomp ) * sr
                  end if
                end do
              end do
            end if
          end if
c preserve parameters of compressed halo(s)
          jcmp = icmp
          do j = 1, ncomp
            icmp = iccmp( j )
            itemp( 1, j ) = imtyp( icmp )
            itemp( 2, j ) = idftyp( icmp )
            stype( j ) = ctype( icmp )
            temp( 2, j ) = cmpmas( icmp )
            temp( 3, j ) = rscale( icmp )
            sdft( j ) = cdft( icmp )
            do i = 1, 3
              temp( 4 + i, j ) = dfcns( i, icmp )
            end do
          end do
c initialize uncompressed halo(s) & DF(s)
          do j = 1, ncomp
            icmp = iccmp( j )
            disc( icmp ) = .false.
            imtyp( icmp ) = ihalo0( j )
            ctype( icmp ) = htype( ihalo0( j ) )
c            cmpmas( icmp ) = temp( 2, j )
c            rscale( icmp ) = temp( 3, j )
            dfcns( 3, icmp ) = cc0( j )
            idftyp( icmp ) = idfn0( j )
            cdft( icmp ) = dftype( idfn0( j ) )
            dfcns( 1, icmp ) = dfm0( j )
            dfcns( 2, icmp ) = dfb0( j )
          end do
c set up uncompressed model
          icmp = iccmp( 1 )
          call inimod
          snglph = .true.
          do j = 1, ncomp
            icmp = iccmp( j )
            r = rtrunc( icmp )
            Emaxo( j ) = Phihal( r )
          end do
          icmp = iccmp( 1 )
          call cutoff( rtrunc( icmp ) )
c set icomp
          icomp = 0
          do i = 1, ncomp
            if( icmp .eq. iccmp( i ) )icomp = i
          end do
          if( icomp .le. 0  )call crash( 'COMPIN', 'nonsense icomp' )
c build or check tables of E_0( J1, J2 ) and df( E, Lz )
          do icomp = 1, ncomp
            icmp = iccmp( icomp )
            x = E0tab( 0.d0, 0.d0 )
            x = max( x, -1.d0 )
            x = distfn( x, 0.d0 )
          end do
c restore parameters of compressed halo(s)
          do j = 1, ncomp
            icmp = iccmp( j )
            imtyp( icmp ) = itemp( 1, j )
            idftyp( icmp ) = itemp( 2, j )
            ctype( icmp ) = stype( j )
            cmpmas( icmp ) = temp( 2, j )
            rscale( icmp ) = temp( 3, j )
            cdft( icmp ) = sdft( j )
            do i = 1, 3
              dfcns( i, icmp ) = temp( 4 + i, j )
            end do
          end do
          cmprssd = .true.
          snglph = .false.
c initialize spline tables
          allocate ( plamda( mradad + 4 ) )
          allocate ( phinc( mradad + 4 ) )
          allocate ( glamda( mradad + 4, ncomp ) )
          allocate ( gmnc( mradad + 4, ncomp ) )
          allocate ( rlamda( mradad + 4, ncomp ) )
          allocate ( rhonc( mradad + 4, ncomp ) )
          r = arad( 1 )
          x = splint2( arad, phin, nradad, r, plamda, phinc, .true. )
          do j = 1, ncomp
            x = splint2( arad, gmn( 1, j ), nradad, r, glamda( 1, j ),
     +                   gmnc( 1, j ), .true. )
            x = splint2( arad, rhon( 1, j ), nradad, r, rlamda( 1, j ),
     +                   rhonc( 1, j ), .true. )
          end do
c reset flag
          iuse = irun
          icmp = jcmp
        end if
        icomp = ncomp
c output part
      else
        if( ( icmp .lt. 1 ) .or. ( icmp .gt. ncmp )
     +                )call crash( 'COMPIN', 'Nonsense value for icmp' )
        ih = imtyp( icmp )
c select file name to use
        if( ih .eq. 27 )then
          fname = 'VKA1.dat'
        else
          fname = 'compress.dat'
        end if
        n = nchars( fname )
c create data file
        open( unit, file = fname( 1:n ), form = 'unformatted',
     +        status = 'new', iostat = i )
        if( i .ne. 0 )then
          if( master )print *, 'compin: .dat file already exists'
          if( gtlogl( 'Do you want to overwrite it' ) )then
            open( unit, file = fname( 1:n ), form = 'unformatted',
     +            status = 'old' )
          else
            call crash( 'COMPIN', '.dat file not overwritten' )
          end if
        end if
c write header records
        write( unit )nradad, ncomp, nrigid, ihalo0, idfn0, iccmp,
     +               ircmp, cc0, dfm0, dfb0, r0core, r1cusp
        do i = 1, nrigid
          j = ircmp( i )
          ldsc( i ) = disc( j )
          jmtyp( i ) = imtyp( j )
        end do
        write( unit )( rtrunc( i ), cmpmas( i ),
     +                 rscale( i ), fmass( i ), i = 1, ncmp ),
     +                          ( ldsc( i ), jmtyp( i ), i = 1, nrigid )
c write out data
        write( unit )( arad( i ), i = 1, nradad ),
     +               ( phin( i ), i = 1, nradad ),
     +               ( ( gmn( i, j ), i = 1, nradad ), j = 1, ncomp ),
     +               ( ( rhon( i, j ), i = 1, nradad ), j = 1, ncomp )
        close( unit )
      end if
      return
      end
