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
c if open of an 'old' file failed, try a different directory
          if( ( istat .ne. 0 ) .and. ( stat .eq. 'old' ) )then
            fname = '../'//run( 4:n )//'/'//run( 1:n )//'.'//inp
            i = lnblnk( fname )
            open( iunit, file = fname( 1:i ),
     +            status = stat, form = frm, iostat = istat )
          end if
        end if
      end if
c report
      if( master )then
        if( istat .ne. 0 )then
          print '( ''Warning: failed to open file '', a, ' //
     +   ' '' on unit'', i3, '' with status = '', a, '', form = '', ' //
     +   ' a, '', access = '', a )', fname( 1:i ), iunit, stat, frm, acc
c        else
c          print '( ''Opened file '', a,  '' on unit'', i3, ' //
c     +          ' '' with status = '', a, '', form = '', a, ' //
c     +      ' '', access = '', a )', fname( 1:i ), iunit, stat, frm, acc
        end if
      end if
      return
      end
