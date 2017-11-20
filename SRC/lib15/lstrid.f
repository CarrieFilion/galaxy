      subroutine lstrid( type, jstep, ilast )
c  Copyright (C) 2015, Jerry Sellwood
c
c To assist the user to select an approriate stride through the .res file
c   when plotting, and to set the last step to be plotted
      implicit none
c
c calling arguments
      character*4 type
      integer ilast, jstep
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
c external
      logical gtlogl
c
c local variables
      integer ifail, ifrst, j
      logical OK
      real am, tms
c
c restart file and read first requested record
      call rewnd
      call nextrec( type, ifail )
      if( ifail .eq. 0 )then
        ifrst = istep
c determine analysis step interval
        call nextrec( type, ifail )
        if( ifail .eq. 0 )then
          am = ts * real( istep - ifrst )
          OK = .false.
          do while ( .not. OK )
            call gtreal(
     +           'Enter time interval for plotting, or 0 for all', tms )
            j = nint( tms / am )
            jstep = j * ( istep - ifrst )
            jstep = max( jstep, 1 )
            if( tms .gt. 0. )then
c check that this is OK
              OK = abs( ts * real( jstep ) - tms ) .lt. .01 * tms
              if( .not. OK )then
                print
     +'( '' Analysis steps every'', f10.2 )', ts * real( istep - ifrst )
                print
     +    '( '' Nearest available stride'', f10.2 )', ts * real( jstep )
                OK = gtlogl( 'Is this acceptable' )
              end if
            else
              OK = .true.
            end if
          end do
          call gtreal(
     +         'Enter time of last moment required, or 0 for all', tms )
          ilast = nint( tms / ts )
          if( ilast .eq. 0 )ilast = 10000000
        else
          print *, 'Only 1 moment available'
          jstep = 1
          ilast = max( istep, 0 )
        end if
c restart the file and re-read first record - do I want this?
        call rewnd
        call nextrec( type, ifail )
      else
        print *, 'No data of requested type ' // type
        jstep = -1
        ilast = -1
      end if
      return
      end
