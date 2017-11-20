      subroutine cmprun( krun )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c Routine to switch to a new .res file for analysis - useful when making
c   comparisons
c
c calling argument
      integer krun
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/anlys.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/comprns.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'inc/supp.f'
c
c local variables
      integer i, j, jcmp, jg, jrun
      logical ldsk
      save ldsk
c
c get list of comparison runs
      if( krun .lt. 0 )then
        nruns = 0
        jrun = irun
        do while ( jrun .gt. 0 )
          nruns = nruns + 1
          if( nruns .gt. mruns )call crash( 'CMPRUN', 'Too many runs' )
          lrun( nruns ) = jrun
          call gtintg(
     +            'Enter new comparison run number, or 0 to end', jrun )
        end do
      else if( krun .le. nruns )then
        if( nruns .gt. 1 )then
          jcmp = icmp
          if( krun .eq. 2 )ldsk = disc( icmp )
c switch to a new file
          close( in )
          irun = lrun( krun )
          jg = jgrid
          call opnfil( in, 'res', 'unformatted', 'old', 'seq', i )
          if( i .ne. 0 )call crash( 'CMPRUN', '.res file not found' )
          rewind in
c read header records
          call hedrc1( in, .true. )
          call codeset
c grid information
          do j = 1, ngrid
            call switch( j )
            call hedrcg( in, .true. )
            if( s3d )then
              i = nr( j )
              s3rad( i ) = rgrid( j )
            end if
            if( p2d .or. p3d )stdpol = ( uoffs .eq. 0. )
            call setgrd( .false. )
          end do
c set up information
          if( rngs )ncmp = ncmp - 1
          do icmp = 1, ncmp
            call hedrcc( in, .true. )
          end do
          icmp = jcmp
c allow for the disk being a different component - special code for Joel!
          if( disc( icmp ) .neqv. ldsk )then
            do i = 1, ncmp
              if( disc( i ) .eqv. ldsk )icmp = i
            end do
          end if
c set number of header records
          nhead = 1 + ncmp
          do j = 1, ngrid
            call switch( j )
            if( .not. ( noslfg .or. dr3d ) )nhead = nhead + 1
          end do
c restore old grid pointer
          if( ngrid .gt. 1 )call switch( jg )
        end if
      else
        print *, 'krun, nruns =', krun, nruns
        call crash( 'CMPRUN', 'calling argument krun > nruns' )
      end if
      return
      end
