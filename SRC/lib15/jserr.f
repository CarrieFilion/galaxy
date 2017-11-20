      subroutine jserr
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c an almost useless error trapping routine - no information is
c   provided as to why it was called or where the call came from!
c
      call jsebuf
      print *,'JSERR called'
      call jsend
      stop
      end
