      integer function nchars( a )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c counts the characters up to the last non-blank one in a string
c   this is redundantly equivalent to the Fortran utility lnblnk
c
c calling argument
      character*(*) a
c
c local variable
      logical blank
c
      nchars = len( a )
      blank = .true.
      do while ( blank )
        if( nchars .gt. 0 )then
          blank = a( nchars:nchars ) .eq. ' '
          if( blank )nchars = nchars - 1
        else
          blank = .false.
        end if
      end do
      return
      end
