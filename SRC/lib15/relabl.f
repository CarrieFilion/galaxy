      subroutine relabl( jst )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to reset grid and plane labels after particles have been moved
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
      include 'inc/model.f'
c
c external
      logical ongrd1
c
c local variables
      integer i, is, j
c
      if( hybrid )then
c sometimes want to force 2nd grid
        if( lhyb2 )then
          do is = 1, jst
            label( is ) = 2
          end do
        else
c determine which grid the particle is on
          do is = 1, jst
            j = igrd( iflag( is ) )
            if( j .eq. 1 )then
c check for particles switching grids
              if( ongrd1( is ) )then
                label( is ) = 1
              else
                label( is ) = 2
              end if
            else
              label( is ) = j
            end if
          end do
        end if
      else if( twogrd )then
c select grid for appropriate population
        do is = 1, jst
          j = igrd( iflag( is ) )
          label( is ) = j
        end do
      else
c only one grid
        do is = 1, jst
          label( is ) = 1
        end do
      end if
c determine new plane number, if it may have changed
      if( nplanes .gt. 1 )then
        j = 0
        do i = 1, ngrid
          if( igrid( i ) .eq. 2 )j = i
        end do
        if( j .eq. 0 )call crash( 'RELABL', 'No c3d active?' )
        do is = 1, jst
          if( ( iz( is ) .lt. nlists ) .and.
     +        ( label( is ) .eq. j ) )then
            ipln( is ) = newc( 3, is ) - xcpred( 3, 1, 1 ) + zm( 1 )
            ipln( is ) = max( ipln( is ), 0 ) + 1
            ipln( is ) = min( ipln( is ), nplanes )
          else
            ipln( is ) = 1
          end if
        end do
      end if
      return
      end
