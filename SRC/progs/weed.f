      program weed
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c copies a runXXX.res file, taking care of sundry errors
c   Unlike merge, there is no option to split the file into the
c   different data types, or to produce a summary listing.   Instead
c   this offers the freedom to select those data types to be copied
c   to the output file and to reduce the number of separate moments
c   of those, if desired.
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
c externals
      character*4 datnam
      logical gtlogl, qcopy
c
c local arrays
      integer ktype
      parameter ( ktype = mtype * mcmp )
      integer jcount( mcmp, mtype ), jstep( mcmp, mtype )
      logical l( mcmp, mtype )
c
c local variables
      integer i, isteps, istp0, itype, j, jcmp, jrun, jstp, k
      integer kstp, nn, nw
      logical cop2, lbad, separ, summ
c
      data l / ktype * .false. /
      data jcount / ktype * 0 /, jstep / ktype * -10 /
c
c set defaults
      call setnag
      master = .true.
      call boilrp( .true. )
c open file and read header
      call header( .false. )
      summ = .false.
c
      print *, 'Enter step intervals between moments to copy for ' //
     +         'both general records and for PNTS, or 0 for all'
      read *, jstp, kstp
      jstp = max( jstp, 1 )
      kstp = max( kstp, 1 )
c write header on output file
      separ = .false.
      do itype = 1, mtype
        do i = 1, mcmp
          l( i, itype ) = .true.
        end do
      end do
      nout = -1
      call opnfil( nout, 'cpy', 'unformatted', 'unknown', 'seq', i )
      call hedrec( nout, .false. )
      isteps = 1
c
      jrun = irun
      read( in, end = 98 )irun, bhol, istep, mr, ma
      print 204, jrun, istep
      istp0 = 0
      if( istep .gt. 0 )then
        print *, 'first step of this file is', istep
        if( gtlogl(
     +         'Do you wish to edit the new file to start from step 0' )
     +                                                    )istp0 = istep
      end if
c rewind and skip forward over header
      call rewnd
c read next header record
   99 read( in, end = 98, err = 91 )irun, bhol, istep, mr, ma
   96 jcmp = 1
      if( irun .le. 10 )then
        jcmp = irun
        if( jcmp .gt. mcmp )call crash( 'WEED', 'jcmp > mcmp' )
      else if( irun .ne. jrun )then
        print *, 'Change of run no from', jrun, ' to', irun
        if( gtlogl( 'Do you want to continue?' ) )go to 96
        stop
      end if
      istep = istep - istp0
      call decode( bhol, itype, nw )
      if( itype .eq. 0 )go to 99
      if( itype .gt. mtype )then
        print *, 'mtype arrays too small', itype, ' needed'
        call crash( 'WEED', 'mtype too small' )
      end if
c 
      if( jstep( jcmp, itype ) .lt. 0 )then
        l( jcmp, itype ) = qcopy( bhol )
      end if
      if( separ )then
        nout = 50 + itype
        if( jstep( jcmp, itype ) .lt. 0 )then
c          l( itype ) = qcopy( bhol )
          if( l( jcmp, itype ) )call hedrec( nout, .false. )
        end if
      end if
      cop2 = l( jcmp, itype )
      if( itype .eq. 1 )then
        if( cop2 )cop2 = mod( istep, kstp ) .eq. 0
      else
        if( cop2 .and. 
     +      ( itype .ne. 26 ) )cop2 = mod( istep, jstp ) .eq. 0
      end if
      cop2 = cop2 .and. ( istep .gt. jstep( jcmp, itype ) )
c      if( .not. cop2 )cop2 = irun .ne. jrun
      if( l( jcmp, itype ) .and. ( .not. cop2 ) )
     +                      print 205, bhol, istep, jstep( jcmp, itype )
c increment size of wres, if necessary
      k = size( wres )
      if( nw .gt. k )then
        if( k .gt. 0 )deallocate ( wres )
        allocate( wres( nw ) )
      end if
c picture data
      if( itype .eq. 1 )then
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
        jstep( jcmp, 1 ) = max( istep, jstep( jcmp, 1 ) )
        if( .not. cop2 )go to 99
        jcount( jcmp, 1 ) = jcount( jcmp, 1 ) + 1
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
c general form of record
      if( cop2 )then
        read( in, end = 98, err = 91 )( wres( i ), i = 1, nw )
        if( ( itype .eq. 11 ) .and. ( ma .eq .ng ) )
     +                                        call dancomp( mr, ma, nw )
        write( nout )irun, bhol, istep, mr, ma
        write( nout )( wres( i ), i = 1, nw )
        jcount( jcmp, itype ) = jcount( jcmp, itype ) + 1
      else
        read( in, end = 98 )
      end if
      jstep( jcmp, itype ) = max( istep, jstep( jcmp, itype ) )
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
        do jcmp = 1, mcmp
          if( jstep( jcmp, itype ) .ge. 0 )print 203,
     +                     jcount( jcmp, itype ), datnam( itype ), jcmp,
     +                                              jstep( jcmp, itype )
        end do
      end do
      close( unit = in )
      stop
  203 format( i5, 2x, a4, ' records copied for compnt', i3,
     +                            ' last complete record for step', i8 )
  204 format( 5x, 'Run no', i6, ' first step no of new file', i8 )
  205 format( 3x, a4, ' record type skipped for step no', i8,
     +        ' previous highest step was', i8 )
  206 format( ' Error reading data record on unit', i4,
     +        ' control record was' )
  207 format( i6, 2x, a4, i8, 2i5 )
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
