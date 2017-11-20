      subroutine jschsz( chsize )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set the size of characters to be used in output strings
c   for historical reasons, the default value of chsize is 0.3
c
c calling argument
      real chsize
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c store new character size (assumed in cm)
      height = chsize
c
c pgplot
      call pgsch( 3. * height )
      return
      end
