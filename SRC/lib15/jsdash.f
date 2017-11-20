      subroutine jsdash( i, j, k, l )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to vary the current line style.  If j=0 the line is continuous
c   otherwise the lengths of line segments are specified by i & k, while
c   the sizes of gaps are specified by j & l
c
c calling arguments
      integer i, j, k, l
c
c common blocks
c
      include 'inc/jscmmn.f'
c
c local arrays
      integer il( 4 )
c
c local variables
      integer m
c
      dash = j .gt. 0
      if( dash )then
        il( 1 ) = i
        il( 2 ) = j
        il( 3 ) = k
        il( 4 ) = l
c store new dash lengths - given in mm
        total = 0
        do m = 1, 4
          il( m ) = max( il( m ), 0 )
          segl( m ) = .1 * real( il( m ) )
          if( mod( m, 2 ) .ne. 0 )then
c ensure dot is as long as it is wide
            if( .not. screen )segl( m ) = max( segl( m ), .02 )
          end if
          total = total + segl( m )
        end do
        partd = 0.
        iseg = 1
      end if
      return
      end
