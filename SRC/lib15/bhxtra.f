      subroutine bhxtra( jst )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c adds local test particles to the Barnes-Hut tree
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local arrays
      real cmax( 3 ), cmin( 3 )
c
c local variables
      integer i, is, jt
      real cbound
c
c work over particles
      do is = 1, jst
c cell boundaries
        do i = 1, 3
          cmin( i ) = bound( 1, i )
          cmax( i ) = bound( 1, i ) + tsize
        end do
c start at stem of tree
        jt = 1
    1   if( itup( jt ) .le. 0 )then
c place in either an empty cell or one with just one other particle
          loc( is ) = -jtmax - 1
          ltree( is ) = jt
        else
c work down chain
          jt = itup( jt )
c find appropriate sub-cell
          do i = 1, 3
            cbound = .5 * ( cmin( i ) + cmax( i ) )
            if( ( oldc( i, is ) - cbound ) .gt. 0. )then
              jt = jt + i
              if( i .eq. 3 )jt = jt + 1
              cmin( i ) = cbound
            else
              cmax( i ) = cbound
            end if
          end do
          go to 1
        end if
      end do
      return
      end
