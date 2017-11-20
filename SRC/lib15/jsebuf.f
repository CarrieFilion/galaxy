      subroutine jsebuf
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to empty the current plotting buffer, to ensure that all
c   lines, text, etc that have been created appear on the device
c   after this call
c
c pgplot
      call pgebuf
      return
      end
