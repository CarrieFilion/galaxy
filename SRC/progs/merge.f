      program merge
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c Copies .res file, correcting sundry errors and optionally creates a
c   printable summary file
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      character*4 datnam
      logical gtlogl
c
c local arrays
      integer ktype
      parameter ( ktype = mtype * mcmp )
      integer jcount( mcmp, mtype ), jstep( mcmp, mtype )
c
c local variables
      character*4 type
      integer i, ip, itype, iw, j, jrun, k, nn, nsum, nw
      logical copy, cop2, lbad, summ
      real c( 3 ), f( 3 ), h( 3 ), ket, pet, p( 3 ), tm, w
      equivalence ( w, iw )
c
      data jcount / ktype * 0 /, jstep / ktype * -10 /
c
c set defaults
      call setnag
      master = .true.
      call boilrp( .true. )
c open file and read header
      call header( .false. )
c
      copy = .true.
c decide about summary file
      summ = gtlogl( 'Do you want to create a summary file?' )
      if( summ )then
        nsum = -1
        call opnfil( nsum, 'sum', 'formatted', 'unknown', 'seq', i )
c        call gtintg( 'Enter nmber of steps between analyses', isteps )
        copy = gtlogl( 'Do you also want to copy the file' )
      end if
c copy header
      if( copy )then
        nout = -1
        call opnfil( nout, 'cpy', 'unformatted', 'unknown', 'seq', i )
        call hedrec( nout, .false. )
      end if
c
      jrun = irun
      read( in, end = 98 )irun, bhol, istep, mr, ma
      print 204, jrun, istep
      if( summ )then
        if( twod )write( nsum, 201 )jrun
        if( threed )write( nsum, 202 )jrun
      end if
c rewind and skip forward over header
      call rewnd
c read next header record
   99 read( in, end = 98, err = 91 )irun, bhol, istep, mr, ma
   96 icmp = 1
      if( irun .le. 10 )then
        icmp = irun
        if( icmp .gt. mcmp )call crash( 'MERGE', 'icmp > mcmp' )
      else if( irun .ne. jrun )then
        print *, 'Change of run no from', jrun, ' to', irun
        if( gtlogl( 'Do you want to continue?' ) )go to 96
        stop
      end if
      call decode( bhol, itype, nw )
      if( itype .eq. 0 )go to 99
      if( itype .gt. mtype )then
        print *, 'mtype arrays too small', itype, ' needed'
        call crash( 'MERGE', 'mtype too small' )
      end if
c
      cop2 = istep .gt. jstep( icmp, itype )
      if( .not. cop2 )print 205, bhol, istep, jstep( icmp, itype )
      cop2 = copy .and. cop2
      write( type, '( a4 )' )bhol
c increment size of wres, if necessary
      k = size( wres )
      if( nw .gt. k )then
        if( k .gt. 1 )deallocate ( wres )
        allocate( wres( nw ) )
      end if
c picture data
      if( type .eq. 'PNTS' )then
c        cop2 = cop2 .and. mod( istep, 200 ) .eq. 0
        k = 0
        nn = mr
        nw = ma
        do while ( nw .gt. 0 )
          j = k + 1
          k = k + nw
          if( cop2 )then
            read( in, end = 98, err = 91 )( wres( i ), i = j, k )
          else
            read( in, end = 98 )
          end if
          nn = nn - nw
          nw = min( nn, nw )
        end do
        jstep( icmp, 1 ) = max( istep, jstep( icmp, 1 ) )
        if( .not. cop2 )go to 99
        jcount( icmp, 1 ) = jcount( icmp, 1 ) + 1
        write( nout )irun, bhol, istep, mr, ma
        k = 0
        nn = mr
        nw = ma
        do while ( nw .gt. 0 )
          j = k + 1
          k = k + nw
          write( nout )( wres( i ), i = j, k )
          nn = nn - nw
          nw = min( nn, nw )
        end do
        go to 99
      end if
c conservation integrals etc.
      if( type .eq. 'INTG' )then
        print 207, irun, bhol, istep, mr, ma
        nw = mr
        read( in, end = 98, err = 91 )( wres( i ), i = 1, nw )
c
        if( summ .and. ( istep .gt. jstep( icmp, itype ) ) )then
          do i = 1, 3
            c( i ) = 0
            p( i ) = 0
            f( i ) = 0
            h( i ) = 0
          end do
          pet = 0
          ket = 0
c sum over active components
          j = 4
          k = 0
          tm = 0
          do ip = 1, ma
            if( uqmass )then
              tm = tm + wres( j + 1 )
              do i = 1, 3
                c( i ) = c( i ) + wres( j + i + 1 ) * wres( j + 1 )
                p( i ) = p( i ) + wres( j + i + 4 ) * wres( j + 1 )
                f( i ) = f( i ) + wres( j + i + 7 )
                h( i ) = h( i ) + wres( j + i + 10 )
              end do
            else
              w = wres( j + 1 )
              k = k + iw
              do i = 1, 3
                c( i ) = c( i ) + wres( j + i + 1 ) * iw
                p( i ) = p( i ) + wres( j + i + 4 ) * iw
                f( i ) = f( i ) + wres( j + i + 7 )
                h( i ) = h( i ) + wres( j + i + 10 )
              end do
            end if
            pet = pet + wres( j + 14 )
            ket = ket + wres( j + 15 )
            j = j + 15
          end do
          w = wres( 4 )
          noff = iw
          do i = 1, 3
c convert CoM position to grid units
            if( uqmass )then
              c( i ) = c( i ) / ( lscale * tm )
            else
              c( i ) = c( i ) / ( lscale * real( k ) )
            end if
c convert drift rate to total momentum in external units
            p( i ) = p( i ) * fmass( 1 ) * gvfac / real( nbod )
          end do
          if( twod )write( nsum, 209 )istep, wres( 2 ), pet, ket,
     +           pet + ket, h( 1 ), p( 1 ), p( 2 ), c( 1 ), c( 2 ),
     +           f( 1 ), f( 2 ), noff
          if( threed )write( nsum, 210 )istep, wres( 2 ), pet, ket,
     +           pet + ket, h, p, c, f, noff
        end if
        if( cop2 )then
          write( nout )irun, bhol, istep, mr, ma
          write( nout )( wres( i ), i = 1, nw )
          jcount( icmp, itype ) = jcount( icmp, itype ) + 1
        end if
        jstep( icmp, itype ) = max( istep, jstep( icmp, itype ) )
        go to 99
      end if
c general form of record
      if( cop2 )then
        read( in, end = 98, err = 91 )( wres( i ), i = 1, nw )
        if( ( type .eq. 'DANL' ) .and. ( ma .eq. ng ) )
     +                                        call dancomp( mr, ma, nw )
c        if( ( type .eq. 'RHOR' ) .and. ( mr .eq. 20000 ) )then
c          call rhorcp( nw )
c        end if
        write( nout )irun, bhol, istep, mr, ma
        write( nout )( wres( i ), i = 1, nw )
        jcount( icmp, itype ) = jcount( icmp, itype ) + 1
      else
        read( in, end = 98 )
      end if
      jstep( icmp, itype ) = max( istep, jstep( icmp, itype ) )
      go to 99
c error during read
   91 print 206, in
      print 207, irun, bhol, istep, mr, ma
      if( .not. gtlogl( 'Do you want to continue?' ) )then
c try reading more records to see whether a good one can be found
        lbad = .true.
        do while ( lbad )
          read( in, err = 91, end = 98 )irun, bhol, istep, mr, ma
          print *, 'Next control record is'
          print 207, irun, bhol, istep, mr, ma
          lbad = .not. gtlogl( 'Is this any good?' )
        end do
        go to 96
      end if
c end of file
   98 print *, 'end of information'
      do itype = 1, mtype
        do icmp = 1, ncmp
          if( jstep( icmp, itype ) .ge. 0 )print 203,
     +jcount( icmp, itype ), datnam( itype ), icmp, jstep( icmp, itype )
          jcount( icmp, itype ) = 0
        end do
      end do
      close( unit = in )
      stop
  201 format( 1h1, 40x, 6hRun no, i5 // 3x, 4hStep, 2x, 8hVirial C, 5x,
     +        5hPot E, 5x, 5hKin E, 5x, 5hTot E, 4x, 5hTot h, 7x,
     +        14hLinear momenta, 6x, 14hCentre of mass, 6x,
     +        14hForce compnnts, 2x, 7hEscapes / )
  202 format( 1h1, 40x, 6hRun no, i5 // 1x, 4hStep, 1x, 8hVirial C, 3x,
     +        5hPot E, 3x, 5hKin E, 3x, 5hTot E, 2x,
     +        22hAngular momentum comps, 2x,
     +        14hLinear momenta, 6x, 14hCentre of mass, 5x,
     +        14hForce compnnts, 3x, 5hEscps / )
  203 format( i5, 2x, a4, i3,
     +        ' records copied; last complete record for step', i8 )
  204 format( 5x, 'Run no', i6, ' first step no of new file', i8 )
  205 format( 3x, a4, ' record type skipped for step no', i8,
     +        ' previous highest step was', i8 )
  206 format( ' Error reading data record on unit', i4,
     +        ' control record was' )
  207 format( i6, 2x, a10, i8, 2i5 )
  209 format( i7, 4f10.4, f10.5, 6f10.3, i7 )
  210 format( i6, 4f8.3, 3f8.3, 3f6.2, 3f7.2, 3f6.1, i7 )
      end

      subroutine rhorcp( nw )
c routine to cut the size of the RHOR data by 5
      use aarrays
      implicit none
c
c calling argument
      integer nw
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      real fac, r( 6 )
      integer i, j, k
c
      k = 0
      j = 0
      r( 1 ) = 0
      fac = 5 * wres( 2 ) * ( 2. * wres( 1 ) )**3
c compute next 5 radii
    1 do i = 1, 5
        r( i + 1 ) = 2 * wres( j + 2 * i - 1 ) - r( i )
      end do
c compute new results
      wres( k + 1 ) = .5 * ( r( 1 ) + r( 6 ) )
      wres( k + 2 ) = fac / ( r( 6 )**3 - r( 1 )**3 )
      r( 1 ) = r( 6 )
      j = j + 10
      k = k + 2
      if( j .lt. nw )go to 1
      print *, j, k
      mr = mr / 5
      ma = ma * 5
      nw = mr
      return
      end

      subroutine opnfil( iunit, inp, frm, stat, acc, istat )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to open a new file with name runXXX.inp with the specified
c   attributes.  The number XXX is constructed from the run number and
c   the file extension type "inp" is hardwired to be 3 characters.
c
c calling arguments
      character inp*3, frm*(*), stat*(*), acc*(*)
      integer istat, iunit
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      integer lnblnk
c
c local variables
      character fname*48, form*6, proc*4, run*8, t*5
      integer i, j, itime, n
      real a
c
      data run( 1:3 ) / 'run' /, form / '( i0 )' /
c
      if( iunit .le. 0 )call new_unit( iunit )
c construct ASCII string for run number
      a = irun
      if( irun .gt. 0 )a = log10( a + .5 )
      n = 4 + int( a )
      if( n .gt. 9 )call crash( 'OPNFIL', 'Run number out of range' )
      write( form( 4:4 ), '( i1 )' )n - 3
      write( run( 4:n ), form )irun
c branch to use time extension to .dmp, pcs or .act files
      if( ( inp .eq. 'dmp' ) .or. ( inp .eq. 'pcs' ) .or.
     +    ( inp .eq. 'act' ) )then
        itime = nint( real( istep ) * ts )
        a = itime
        if( itime .gt. 0 )a = log10( a + .5 )
        i = 1 + int( a )
        if( i .gt. 5 )call crash( 'OPNFIL', 'time out of range' )
        write( form( 4:4 ), '( i1 )' )i
        write( t( 1:i ), form )itime
c encode processor id in file name
        if( ( inp .eq. 'dmp' ) .and. parallel )then
          a = myid
          if( myid .gt. 0 )a = log10( a + .5 )
          j = 1 + int( a )
          if( ( j .le. 0 ) .or. ( j .gt. 4 ) )call crash( 'OPNFIL',
     +                                     'Processor id out of range' )
          write( form( 4:4 ), '( i1 )' )j
          write( proc( 1:j ), form )myid
          fname = run( 1:n )//'.'//proc( 1:j )//inp//t( 1:i )
          i = lnblnk( fname )
          open( iunit, file = fname( 1:i ),
     +          status = stat, form = frm, iostat = istat )
        else
c single processor option
          fname = run( 1:n )//'.'//inp//t( 1:i )
          i = lnblnk( fname )
          open( iunit, file = fname( 1:i ),
     +          status = stat, form = frm, iostat = istat )
        end if
      else
c other files
        fname = run( 1:n )//'.'//inp
        i = lnblnk( fname )
        if( acc .eq. 'append' )then
          open( iunit, file = fname( 1:i ), status = stat,
     +          form = frm, access = acc, iostat = istat )
        else
          open( iunit, file = fname( 1:i ),
     +          status = stat, form = frm, iostat = istat )
c
c force use of same directory for safety reasons else I could destroy
c   a valuable .res file with a subsequent
c
c       mv fort.3 runXXXX.res
c
c if open of an 'old' file failed, try a different directory
c          if( ( istat .ne. 0 ) .and. ( stat .eq. 'old' ) )then
c            fname = '../'//run( 4:n )//'/'//run( 1:n )//'.'//inp 
c            i = lnblnk( fname )
c            open( iunit, file = fname( 1:i ),
c     +            status = stat, form = frm, iostat = istat )
c          end if
        end if
      end if
c report
      if( master )then
        if( istat .ne. 0 )then
          print '( ''Failed to open file '', a,  '' on unit'', i3, ' //
     +          ' '' with status = '', a, '', form = '', a, ' //
     +      ' '', access = '', a )', fname( 1:i ), iunit, stat, frm, acc
c        else
c          print '( ''Opened file '', a,  '' on unit'', i3, ' //
c     +          ' '' with status = '', a, '', form = '', a, ' //
c     +      ' '', access = '', a )', fname( 1:i ), iunit, stat, frm, acc
        end if
      end if
      return
      end

      real*8 function dfwght( E, Lz )
c function to weight DF when unequal particle masses are requested
c  value of this function is the inverse of the desired particle weights
      implicit none
c
c calling arguments
      real*8 E, Lz
c
c dummy version for analysis only
      dfwght = 1.
      return
      end
