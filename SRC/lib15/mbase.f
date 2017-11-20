      subroutine mbase( m, nmodes, sumsq, search, found )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c part of the mode fitting software
c it maintains a direct access file containing results from fits
c   it will add new records to the file and check whether a
c   given set of data have previously been fitted with precisely
c   the same parameters
c
c calling arguments
      integer m, nmodes
      logical search, found
      real sumsq
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      integer mrun
      parameter ( mrun = 100 )
      common / savfrq / imb, nrec, nrun, run( mrun )
      integer imb, nrec, nrun, run
c
c local variables
      character*4 bstr( 5 )
      integer i, irec, jjp, jjt, jrun, kkp, kkt, modes, n
      real ahol, c, frq( 2, mmodes ), sum, type( 5 )
c
      data imb / 0 /
      data bstr / 'LGSP', 'DANL', 'SFPC', 'SPHB', 'ZANL' /
c
c skip if file is already open
      if( imb .ne. irun )then
c try to open an old file
        i = 40 + 8 * mmodes
        open( unit = 13, err = 2, file = '../mbase',
     +        access = 'direct', recl = i, status = 'old' )
c find last record
        nrec = 0
        nrun = 0
        jrun = 1
        do while ( jrun .gt. 0 )
          nrec = nrec + 1
          read( 13, rec = nrec )jrun, ahol, n, jjt, kkt, jjp, kkp,
     +                          modes, c, sum, frq
          if( irun .eq. jrun )then
            nrun = nrun + 1
            if( nrun .gt. mrun )then
              print *, 'run array size too small - error in MBASE'
              stop
            end if
            run( nrun ) = nrec
          end if
        end do
        nrec = nrec - 1
        print *, nrec, ' records found on old data file', nrun,
     +          ' for this run'
        imb = irun
        return
c open a new file
    2   print *, 'No old file found'
        call gtintg( 'Enter 0 to open a new file, 1 to stop', i )
        if( i .eq. 1 )stop
        print *, 'Opening new file'
        i = 40 + 8 * mmodes
        open( unit = 13, file = '../mbase', access = 'direct',
     +        recl = i, status = 'new' )
c flag last record
        i = -1
        write( 13, rec = 1 )i
        nrec = 0
        imb = irun
        return
      end if
c
c search for existing fit
c
      found = .false.
      if( nrun .gt. 0 )then
        do i = 1, 5
          read( bstr( i ), '( a4 )' )type( i )
        end do
        do i = 1, nrun
          irec = run( i )
          read( 13, rec = irec )jrun, ahol, n, jjt, kkt, jjp, kkp,
     +                          modes, c
          if( jrun .ne. irun )then
            print *, 'logical error in MBASE'
            print *, irun, jrun, i, irec, nrun
            stop
          end if
          if( ( lgsp .and. ( ahol .eq. type( 1 ) ) ) .or.
     +        ( danl .and. ( ahol .eq. type( 2 ) ) ) .or.
     +        ( sfpc .and. ( ahol .eq. type( 3 ) ) ) .or.
     +        ( sphb .and. ( ahol .eq. type( 4 ) ) ) .or.
     +        ( zanl .and. ( ahol .eq. type( 5 ) ) ) )then
            if( ( n .eq. m ) .and.
     +          ( modes .eq. nmodes ) .and.
     +          ( jjt .eq. jt ) .and.
     +          ( kkt .eq. kt ) .and.
     +          ( jjp .eq. jp ) .and.
     +          ( kkp .eq. kp ) .and.
     +          ( c .eq. apodc ) )then
c success!
              found = .true.
              if( search )then
                print *, 'Old solution found'
                read( 13, rec = irec )irun, ahol, m, jt, kt, jjp, kkp,
     +                                nmodes, apodc, sumsq, freq
                return
              else
                go to 1
              end if
            end if
          end if
        end do
      end if
c no matching data found
      if( search )return
      irec = nrec + 1
c save new fit if requested
    1 if( found )then
        print *, 'Over-writing old solution'
      else
        i = 0
        if( lgsp )i = 1
        if( danl )i = 2
        if( sfpc )i = 3
        if( sphb )i = 4
        if( zanl )i = 5
        if( i .eq. 0 )call crash( 'MBASE', 'Unrecognized data type' )
        read( bstr( i ), '( a4 )' )ahol
      end if
      write( 13, rec = irec )irun, ahol, m, jt, kt, jp, kp, nmodes,
     +                       apodc, sumsq, freq
      if( found )return
c advance last record flag
      nrun = nrun + 1
      if( nrun .gt. mrun )then
        print *, 'run array size too small - error in MBASE'
        stop
      end if
      run( nrun ) = irec
      i = -1
      nrec = nrec + 1
      write( 13, rec = nrec + 1 )i
      return
      end
