      subroutine nextrec( type, ifail )
c  Copyright (C) 2015, Jerry Sellwood
c
c General purpose routine in the analysis software to read the next
c   available record from the .res file.  The record type is specified
c   by a 4 character string, which is interpreted in subroutine DECODE.
c   The record lengths of the various types are also calculated by DECODE.
c
c If the different data types have been split into separate files (log
c   variable 'separt' set to .true.) then the unit number of the requested
c   type is related to itype, and the very next record is read from that
c   file.  Otherwise, the routine reads record headers and skips data
c   records in the main runXXX.res file until the key word in the record
c   header matches the requested type.
c
c The data from the requested record are read into the array wres
      use aarrays
      implicit none
c
c calling arguments
      character*4 type
      integer ifail
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c local variables
      character*4 b
      integer i, ip, itype, iunit, j, k, la, lr, nw
      logical copy
      save ip
c
c allocate tme array on first call
      if( .not. allocated( tme ) )then
        ntm = mtm
        allocate ( tme( ntm ) )
      end if
c find unit number for file
      if( separt )then
        call decode( type, itype, nw )
        iunit = 50 + itype
        if( .not. lunit( itype ) )then
          do i = 1, nhead
            read( iunit )
          end do
        end if
        lunit( itype ) = .true.
      else
        iunit = in
      end if
c preserve last record parameters
      lr = mr
      la = ma
c search for next record of requested type
      ifail = 0
      copy = .false.
      do while ( .not. copy )
c get next record header
        read( iunit, end = 98, err = 98 )ip, bhol, istep, mr, ma
c convert Hollerith string to type character
        write( b, '( a4 )' )bhol
        if( byswap )call swap( b )
        call decode( b, itype, nw )
c an extra file header is flagged as type 0
        if( itype .eq. 0 )then
          call crash( 'NEXTREC', 'Extra header record found' )
        else
          copy = ( type .eq. b )
c allow for ip flagging a grid number
          if( copy .and. lgtyp( itype ) .and. ( ip .lt. mcmp ) )then
            k = -1
            do i = 1, ncmp
              if( igrd( i ) .eq. ip )k = i
            end do
            if( k .le. 0 )call crash( 'NEXTREC',
     +           'Could not match the grid data to any mass component' )
            ip = k
          end if
c check whether ip matches the desired component
          if( ip .lt. mcmp )then
            copy = copy .and. ( ip .eq. icmp )
          else
            ip = 1
          end if
c skip if this record is at the same or an earlier time than one already read
          if( copy )copy = copy .and. ( lstp( ip, itype ) .lt. istep )
          lstp( ip, itype ) = max( istep, lstp( ip, itype ) )
          time = real( istep ) * ts
c increment size of wres, if necessary
          k = size( wres )
          if( nw .gt. k )then
            if( allocated( wres ) )deallocate ( wres )
            allocate( wres( nw ) )
          end if
c read picture data
          if( itype .eq. 1 )then
            k = 0
            do while ( k .lt. mr )
              j = k + 1
              k = k + ma
              k = min( k, mr )
              if( copy )then
                read( iunit, end = 98, err = 98 )( wres( i ), i = j, k )
              else
                read( iunit, end = 98, err = 98 )
              end if
            end do
          else
c check space and read in data
            if( copy )then
              read( iunit, end = 98, err = 98 )( wres( i ), i = 1, nw )
              if( ( itype .eq. 11 ) .and. ( ma .eq. ng )
     +                                       )call dancomp( mr, ma, nw )
            else
              read( iunit, end = 98, err = 98 )
            end if
          end if
        end if
      end do
      return
c end of file
   98 print *, ' eof'
c restore record parameters for last successful read
      mr = lr
      ma = la
      call decode( type, itype, nw )
      read( type, '( a4 )' )bhol
      ip = 1
      do i = 1, mcmp
        if( lstp( i, itype ) .gt. 0 )ip = i
      end do
      istep = lstp( ip, itype )
      time = real( istep ) * ts
      ifail = 1
      return
      end
