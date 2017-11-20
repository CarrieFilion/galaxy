      subroutine mascmb
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to sum masses in different zones on the grid prior to force
c   determination
c Parallel version that combines masses from different nodes, when needed,
c   and also works in single processor mode
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variables
      integer i, ig, k
c
c sum mass arrays from all processors
      if( parallel )call massum
c
c loop over methods
      do ig = 1, ngrid
        call switch( ig )
        if( dr3d )then
c copy ptcl coordinates from new to old positions
          do k = 1, ndrct
            do i = 1, 3
              drpt( i + 1, k ) = drpt( i + 4, k )
            end do
          end do
        end if
      end do
c combine masses from other zones
      call zoncmb
c make predicted centroid position of zone 1 the actual centroid pos
      if( centrd )then
        do ig = 1, ngrid
          do i = 1, 3
            xcen( i, ig ) = xcpred( i, 1, ig )
          end do
        end do
      end if
c sweep over s3d grid rings to tabulate interior and exterior masses
      call s3dswp
      return
      end
