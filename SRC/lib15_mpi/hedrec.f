      subroutine hedrec( iunit, inp )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to write or read the header of .dmp or .res files
c this information is supposed to be enough to permit a run to be
c   restarted and all the variables required for the analysis software
c
c calling arguments - the unit for the file and input/output switch
      integer iunit
      logical inp
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      character*4 dftype, htype
c
c local array
      integer ig( 10 )
c
c local variables
      integer i, j, ndim
      logical l
c
      if( master )rewind iunit
c
c first header record
c
      if( inp )ndim = ndimen
      call hedrc1( iunit, inp )
c set logicals
      if( inp )then
        call codeset
        stdpol = .false.
        if( ngrid .gt. 2 )then
          if( master )print *, ngrid, ( igrid( i ), i = 1, 3 )
          call crash( 'HEDREC', 'More than 2 grids' )
        end if
c check that number of dimensions has not changed
        if( ndim .ne. ndimen )then
          if( master .and. ( ndim .gt. 0 ) )print *,
     +            'ndimen was reset in hedrec from', ndim, ' to', ndimen
        end if
c set some values
        gvfac = 1. / ( lscale * ts )
        nsect = max( nsect, 1 )
        epsiln = softl / lscale
c tsoft could be undefined in old headers, default is 1
        if( tsoft .le. 0 )tsoft = 1
      end if
c
c grid information
c
      do j = 1, ngrid
        call switch( j )
        call hedrcg( iunit, inp )
      end do
      call switch( 0 )
c process grid information
      if( inp )then
c process grid information
        do j = 1, ngrid
          call switch( j )
          if( s3d )then
            i = nr( j )
            s3rad( i ) = rgrid( j )
          end if
          if( p2d .or. p3d )stdpol = ( uoffs .eq. 0. )
          call setgrd( .true. )
        end do
        call switch( 0 )
      end if
c
c component information
c
      do icmp = 1, ncmp
        call hedrcc( iunit, inp )
      end do
      if( inp )then
c process this information
        l = .false.
        rmax = 0
        if( rngs )ncmp = ncmp - 1
        do i = 1, ncmp
          testp( i ) = cmpmas( i ) .eq. 0.
          l = l .or. ( ctype( i ) .eq. 'ADIA' )
          rmax = max( sngl( rmax ), rtrunc( i ) )
          cdft( i ) = dftype( idftyp( i ) )
        end do
c initialize variables in / model /
c read data on the compressed halo
        if( l )then
          do i = 1, ncmp
            if( ctype( i ) .eq. 'ADIA' )icmp = i
          end do
          call compin( .true. )
        end if
        call inimod
        do i = 1, ncmp
          rigidp( i ) = nsp( i ) .eq. 0
        end do
c set truncation and angular momentum tapers
        call cutoff( sngl( rmax ) )
c output some information
        if( master )then
          write( no, '( i10, '' particles in'', i3, '' populations'' )'
     +               )nbod, ncmp
          do i = 1, ncmp
            write( no, '( ''Pop'', i3, '' has mass'' f6.2,' //
     +               ' '' and scale length'', f6.2, '' grid units'' )' )
     +                                       i, cmpmas( i ), rscale( i )
          end do
          write( no, '( 11x, ''Time unit'', f10.3 )' )ts
          if( nsect .gt. 1 )write( no, '( '' Mass distribution is'',' //
     +                                ' i4,''-fold symmetric'' )' )nsect
        end if
c external perturber
        if( pertbn .and. master )then
          if( exsphr .or. exgnrc )then
            if( exgnrc )then
              write( no, * )'Generic perturber'
              write( no, '( '' Perturber mass'', 2f10.4 )' )mptbr
            else
              write( no, * )'Spherical perturber ' // htype( isphp )
              write( no,
     +      '( '' Perturber or central mass and core radius'', 2f10.4 )'
     +         )mptbr, eptbr
            end if
          else if( extbar )then
            write( no,
     +'( '' Perturbing bar mass and major axis'', 2f10.4 )' )mptbr, abar
          else if( exspir )then
            write( no, * )'Spiral perturber'
            write( no,
     +     '( '' Pattern speed and amplitude'', 2f10.4 )' )omegsp, epspi
          else
            call crash( 'HEDREC', 'Unknown perturber 2' )
          end if
        end if
        if( master )then
c print out Fourier filter only if it is set up
          l = .false.
          do i = 1, ng
            l = l .or. lg( i )
          end do
          if( l )then
            write( no, * )'Active sectoral harmonics are:'
            j = 0
            do i = 1, ng
              if( .not. lg( i ) )then
                j = j + 1
                if( j .gt. 10 )then
                  write( no, '( 10i8 )' )ig
                  j = 1
                end if
                ig( j ) = i - 1
              end if
            end do
            write( no, '( 10i8 )' )( ig( i ), i = 1, j )
          end if
        end if
      end if
      return
      end

      subroutine hedrc1( iunit, linp )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c routine to assemble or distribute the information for the first header record
c
c calling arguments
      integer iunit
      logical linp
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
      include 'mpif.h'
c
c local array
      integer idim
      parameter ( idim = 52 )
      integer idm( idim )
      logical ldm( idim )
      real dum( idim )
      equivalence ( dum, idm ), ( dum, ldm )
c
c local variable
      integer i, j
c
      if( linp )then
        if( master )then
          read( iunit )dum
c flag version change
          if( dum( 1 ) .lt. cvers )then
            print '( / 10x, 30(''*'') // 10x, ''Warning: version '' ' //
     +            ' ''of code used when data file was created'', ' //
     +            ' f6.2 / 37x, ''version in this executable is'', ' //
     +            ' f6.2 // 10x, 30(''*'') )', dum( 1 ), cvers
          end if
        end if
        if( parallel )call mpi_bcast( idm, idim, mpi_integer, 0,
     +                                               mpi_comm_world, i )
c distribute information
        cvers = dum(  1 )
        ngrid = idm(  2 )
        do i = 1, 3
          igrid( i ) = idm(  2 + i )
        end do
        lscale = dum(  6 )
        ts     = dum(  7 )
        nlists = idm(  8 )
        ng     = idm(  9 )
        nbod   = idm( 10 )
        ncmp   = idm( 11 )
        pmass  = dum( 12 )
        softl  = dum( 13 )
        fixrad = ldm( 14 )
        nsect  = idm( 15 )
        jmass  = idm( 16 )
        hybrid = ldm( 17 )
        twogrd = ldm( 18 )
        heavies= ldm( 19 )
        kdrct  = idm( 20 )
        tsoft  = idm( 21 )
        ndimen = idm( 22 )
        alp2   = dum( 23 )
        drfrqs = dum( 24 )
c entries 25 thru 33 are unused
        cbar   = dum( 34 )
        alpha2 = dum( 35 )
        beta2  = dum( 36 )
        j      = idm( 38 )
        mptbr  = dum( 39 )
c interpret perturber information
        pertbn = j .ne. 0
        if( pertbn )then
          exsphr = j .lt. 0
          if( exsphr )isphp = -j
          extbar = j .eq. 2
          exspir = j .eq. 4
          exgnrc = j .eq. 5
          if( j .gt. 5 )call crash( 'HEDRC1', 'Unknown perturber' )
          if( exsphr )then
            eptbr   = dum( 40 )
          else if( exspir )then
            epspi   = dum( 40 )
            omegsp  = dum( 37 )
          else if( extbar )then
            abar    = dum( 40 )
            bbar   = dum( 37 )
            eptbr = bbar / abar
          end if
        end if
        centrd = ldm( 41 )
        offfor = ldm( 42 )
        uqmass = ldm( 43 )
        lfrv   = ldm( 44 )
        offanal= ldm( 45 )
        selfe  = dum( 46 )
        drfac  = dum( 47 )
        Lzfac  = dum( 48 )
        nprop  = idm( 49 )
        rbess  = dum( 50 )
        hzfac  = dum( 51 )
        zshift = dum( 52 )
      else
c assemble information
        dum(  1 ) = cvers
        idm(  2 ) = ngrid
        do i = 1, 3
          idm( i + 2 ) = igrid( i )
        end do
        dum(  6 ) = lscale
        dum(  7 ) = ts
        idm(  8 ) = nlists
        idm(  9 ) = ng
        idm( 10 ) = nbod
        idm( 11 ) = ncmp
        dum( 12 ) = pmass
        dum( 13 ) = softl
        ldm( 14 ) = fixrad
        idm( 15 ) = nsect
        idm( 16 ) = jmass
        ldm( 17 ) = hybrid
        ldm( 18 ) = twogrd
        ldm( 19 ) = heavies
        idm( 20 ) = kdrct
        idm( 21 ) = tsoft
        idm( 22 ) = ndimen
        dum( 23 ) = alp2
        dum( 24 ) = drfrqs
c entries 25 thru 33 are unused
        dum( 34 ) = cbar
        dum( 35 ) = alpha2
        dum( 36 ) = beta2
        dum( 39 ) = mptbr
c encode perturber information
        j = 0
        if( pertbn )then
          if( exsphr .or. exgnrc )then
            if( exsphr )j = -isphp
            if( exgnrc )j = 5
            dum( 40 ) = eptbr
          else if( extbar )then
            j = 2
            dum( 40 ) = abar
            dum( 37 ) = bbar
          else if( exspir )then
            j = 4
            dum( 40 ) = epspi
            dum( 37 ) = omegsp
          else
            call crash( 'HEDRC1', 'Unknown perturber' )
          end if
        end if
        idm( 38 ) = j
        ldm( 41 ) = centrd
        ldm( 42 ) = offfor
        ldm( 43 ) = uqmass
        ldm( 44 ) = lfrv
        ldm( 45 ) = offanal
        dum( 46 ) = selfe
        dum( 47 ) = drfac
        dum( 48 ) = Lzfac
        idm( 49 ) = nprop
        dum( 50 ) = rbess
        dum( 51 ) = hzfac
        dum( 52 ) = zshift
        if( master )write( iunit )dum
      end if
      return
      end

      subroutine hedrcg( iunit, linp )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c routine to assemble or distribute the information for the first header record
c
c calling arguments
      integer iunit
      logical linp
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'mpif.h'
c
c local array
      integer idim
      parameter ( idim = 80 )
      integer idm( idim )
      logical ldm( idim )
      real dum( idim )
      equivalence ( dum, idm ), ( dum, ldm )
c
c local variables
      integer i, j
c
      if( p2d .or. p3d .or. p3a .or. s3d )then
        j = ng + 9
        if( j .gt. idim )call crash( 'HEDRCG', 'idim too small' )
        if( linp )then
          if( master )read( iunit )( dum( i ), i = 1, j )
          if( parallel )call mpi_bcast( dum, j, mpi_real, 0,
     +                                               mpi_comm_world, i )
c distribute information
          nr( jgrid )
     +           = idm(  1 )
          na     = idm(  2 )
          ng     = idm(  3 )
          rgrid( jgrid )
     +           = dum(  4 )
          gamma  = dum(  5 )
          uoffs  = dum(  6 )
          ngz    = idm(  7 )
          dzg    = dum(  8 )
          s3lmax = idm(  9 )
          do i = 1, ng
            lg( i ) = ldm( 9 + i )
          end do
        else
c assemble information
          idm(  1 ) = nr( jgrid )
          idm(  2 ) = na
          idm(  3 ) = ng
          dum(  4 ) = rgrid( jgrid )
          dum(  5 ) = gamma
          dum(  6 ) = uoffs
          idm(  7 ) = ngz
          dum(  8 ) = dzg
          idm(  9 ) = s3lmax
          do i = 1, ng
            ldm( 9 + i ) = lg( i )
          end do
          j = ng + 9
          if( master )write( iunit )( dum( i ), i = 1, j )
        end if
      else if( sf2d .or. sf3d )then
        if( linp )then
          if( master )then
c clumsy: lastf is not yet defined!
            read( iunit )( dum( i ), i = 1, 5 )
            lastf = idm( 5 )
            backspace( iunit )
          end if
          if( parallel )call mpi_bcast( lastf, 1, mpi_integer, 0,
     +                                               mpi_comm_world, i )
        end if
        j = 2 * lastf + 7
        if( j .gt. idim )call crash( 'HEDRCG', 'idim too small' )
        if( linp )then
          if( master )read( iunit )( dum( i ), i = 1, j )
          if( parallel )call mpi_bcast( dum, j, mpi_real, 0,
     +                                               mpi_comm_world, i )
c distribute information
          write( basset, '( a4 )' )dum( 1 )
          basis    = idm( 2 )
          minmaxr  = dum( 3 )
          maxmaxr  = dum( 4 )
          lastf    = idm( 5 )
          ngz      = idm( 6 )
          dzg      = dum( 7 )
          j = 7
          do i = 1, lastf
            nsel( i ) = idm( j + 1 )
            msel( i ) = idm( j + 2 )
            j = j + 2
          end do
        else
c assemble information
          read( basset, '( a4 )' )dum( 1 )
          idm( 2 ) = basis
          dum( 3 ) = minmaxr
          dum( 4 ) = maxmaxr
          idm( 5 ) = lastf
          idm( 6 ) = ngz
          dum( 7 ) = dzg
          j = 7
          do i = 1, lastf
            idm( j + 1 ) = nsel( i )
            idm( j + 2 ) = msel( i )
            j = j + 2
          end do
          if( master )write( iunit )( dum( i ), i = 1, j )
        end if
      else if( c2d .or. c3d )then
        if( linp )then
          if( master )read( iunit )( dum( i ), i = 1, 6 )
          if( parallel )call mpi_bcast( dum, 6, mpi_real, 0,
     +                                               mpi_comm_world, i )
c distribute information
          ngx = idm( 1 )
          ngy = idm( 2 )
          ngz
     +         = idm( 3 )
          do i = 1, 3
            dh( i ) = dum( 3 + i )
          end do
        else
c assemble information
          idm( 1 ) = ngx
          idm( 2 ) = ngy
          idm( 3 ) = ngz
          do i = 1, 3
            dum( 3 + i ) = dh( i )
          end do
          if( master )write( iunit )( dum( i ), i = 1, 6 )
        end if
      else if( scf )then
c watch out - how is ng defined?
        if( ng .le. 0 )call crash( 'HEDRCG', 'ng is undefined' )
        j = ng + 6
        if( j .gt. idim )call crash( 'HEDRCG', 'idim too small' )
        if( linp )then
          if( master )read( iunit )( dum( i ), i = 1, j )
          if( parallel )call mpi_bcast( dum, j, mpi_real, 0,
     +                                               mpi_comm_world, i )
c distribute information
          write( basset, '( a4 )' )dum( 1 )
          nbas   = idm( 2 )
          s3lmax = idm( 3 )
          lmin   = idm( 4 )
          lskip  = idm( 5 )
          rgrid( jgrid )
     +           = dum( 6 )
          do i = 1, ng
            lg( i ) = ldm( i + 6 )
          end do
        else
c assemble information
          read( basset, '( a4 )' )dum( 1 )
          idm( 2 ) = nbas
          idm( 3 ) = s3lmax
          idm( 4 ) = lmin
          idm( 5 ) = lskip
          dum( 6 ) = rgrid( jgrid )
          do i = 1, ng
            ldm( i + 6 ) = lg( i )
          end do
          j = ng + 6
          if( master )write( iunit )( dum( i ), i = 1, j )
        end if
      else if( bht )then
        if( linp )then
          if( master )read( iunit )dum( 1 )
          if( parallel )call mpi_bcast( dum, 1, mpi_real, 0,
     +                                               mpi_comm_world, i )
c distribute information
          bhtol = dum( 1 )
        else
c assemble information
          dum( 1 ) = bhtol
          if( master )write( iunit )dum( 1 )
        end if
      else
        call crash( 'HEDRCG', 'Unrecognised metod' )
      end if
      return
      end

      subroutine hedrcc( iunit, linp )
c  Copyright (C) 2016, Jerry Sellwood
      implicit none
c routine to assemble or distribute the information for the first header record
c
c calling arguments
      integer iunit
      logical linp
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      include 'mpif.h'
c
c local array
      integer idim
      parameter ( idim = 38 )
      integer idm( idim )
      logical ldm( idim )
      real dum( idim )
      equivalence ( dum, idm ), ( dum, ldm )
c
c local variables
      integer i, j
c
      if( linp )then
        if( master )read( iunit )dum
        if( parallel )call mpi_bcast( dum, idim, mpi_real, 0,
     +                                               mpi_comm_world, i )
c distribute information
        nsp( icmp )    = idm(  1 )
        cmpmas( icmp ) = dum(  2 )
        rscale( icmp ) = dum(  3 )
        fmass( icmp )  = dum(  4 )
        rtrunc( icmp ) = dum(  5 )
        disc( icmp )   = ldm(  6 )
        write( ctype( icmp ), '( a4 )' )dum( 7 )
        quiet( icmp )  = ldm(  8 ) 
        npr( icmp )    = idm(  9 )
        smr( icmp )    = ldm( 10 )
        dist( icmp )   = ldm( 11 )
        iztyp( icmp )  = idm( 12 )
        z0init( icmp ) = dum( 13 )
        imtyp( icmp )  = idm( 14 )
        impar( icmp )  = idm( 15 )
        idftyp( icmp ) = idm( 16 )
        do i = 1, 3
          jdfcs( i, icmp ) = idm( i + 16 )
          dfcns( i, icmp ) = dum( i + 19 )
        end do
        igrd( icmp )   = idm( 23 )
        do i = 1, 6
          comi( i, icmp ) = dum( i + 23 )
        end do
        Lztmn( icmp )  = dum( 30 )
        Lztmx( icmp )  = dum( 31 )
        Lzcrt( icmp )  = dum( 32 )
        Lzrmn( icmp )  = dum( 33 )
        Lzmno( icmp )  = dum( 34 )
        znorm( icmp )  = dum( 35 )
        DFnorm( icmp ) = dum( 36 )
        cmppar( 1, icmp )
     +                 = dum( 37 )
c                        dum( 38 ) unused
      else
c assemble information
        idm(  1 ) = nsp( icmp )
        dum(  2 ) = cmpmas( icmp )
        dum(  3 ) = rscale( icmp )
        dum(  4 ) = fmass( icmp )
        dum(  5 ) = rtrunc( icmp )
        ldm(  6 ) = disc( icmp )
        read( ctype( icmp ), '( a4 )' )dum( 7 )
        ldm(  8 ) = quiet( icmp )
        idm(  9 ) = npr( icmp )
        ldm( 10 ) = smr( icmp )
        ldm( 11 ) = dist( icmp )
        idm( 12 ) = iztyp( icmp )
        dum( 13 ) = z0init( icmp )
        idm( 14 ) = imtyp( icmp )
        idm( 15 ) = impar( icmp )
        idm( 16 ) = idftyp( icmp )
        do i = 1, 3
          idm( i + 16 ) = jdfcs( i, icmp )
          dum( i + 19 ) = dfcns( i, icmp )
        end do
        idm( 23 ) = igrd( icmp )
        do i = 1, 6
          dum( i + 23 ) = comi( i, icmp )
        end do
        dum( 30 ) = Lztmn( icmp )
        dum( 31 ) = Lztmx( icmp )
        dum( 32 ) = Lzcrt( icmp )
        dum( 33 ) = Lzrmn( icmp )
        dum( 34 ) = Lzmno( icmp )
        dum( 35 ) = znorm( icmp )
        dum( 36 ) = DFnorm( icmp )
        dum( 37 ) = cmppar( 1, icmp )
c       dum( 38 ) unused
        if( master )write( iunit )dum
      end if
      return
      end
