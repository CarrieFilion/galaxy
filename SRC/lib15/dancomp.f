      subroutine dancomp( n, m, nw )
c  Copyright (C) 2015, Jerry Sellwood
      use aarrays
      implicit none
c a little-used routine to compact the information in the sectoral harmomic
c   expansion of the density by selecting only those values that are non-zero
c
c calling arguments
      integer m, n, nw
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local array
      real t( 5000 )
c
c local variables
      integer i, j, k1, k2, l1, l2, mnew
c
      if( nsect .eq. 1 )return
      mnew = ( m - 1 ) / nsect + 1
c      print *, 'dancomp called,  m_new  = ', mnew, ' m_old  = ', m
c copy danl array omitting unwanted values
      l1 = 0
      l2 = n * m
      k1 = 0
      k2 = n * mnew
      do i = 1, n
        do j = 1, m
          l1 = l1 + 1
          l2 = l2 + 1
          if( mod( j - 1, nsect ) .eq. 0 )then
            k1 = k1 + 1
            k2 = k2 + 1
            t( k1 ) = wres( l1 )
            t( k2 ) = wres( l2 )
          end if
        end do
      end do
c overwrite old data
      m = mnew
      nw = 2 * n * m
      do i = 1, nw
        wres( i ) = t( i )
      end do
      return
      end
