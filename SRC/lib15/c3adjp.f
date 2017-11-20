      subroutine c3adjp
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to massage potential on a 3-D Cartesian grid to minimize
c   force anisotropies - not used at present
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
c local array
      real adj( 3 )
c
c local variables
      integer i, j, k, n
      real w
c
c      adj( 1 ) = 0.336505 - 0.451126
c      adj( 2 ) = 0.303310 - 0.451126
c      adj( 3 ) = 0.278240 - 0.352670
c      data adj / -.114621, -.147816, -.074430 /
      adj( 1 ) = 2.83269

      if( .not. c3d )call crash( 'C3ADJP', 'Unrecognised code' )
c assume no mass on boundary planes
      n = ngxy - ngx - 1
      do k = 2, ngz - 1
        n = n + 2 * ngx
        do j = 2, ngy - 1
          n = n + 2
          do i = 2, ngx - 1
            n = n + 1
            w = grdmss( n, 1 )
            if( w .ne. 0. )then
c central value
              grdfld( n, 4 ) = grdfld( n, 4 ) + w * adj( 1 )
c 6 face centres
c              grdfld( n +    1, 4 ) =
c     +                              grdfld( n +    1, 4 ) + w * adj( 2 )
c              grdfld( n -    1, 4 ) =
c     +                              grdfld( n -    1, 4 ) + w * adj( 2 )
c              grdfld( n +  ngx, 4 ) =
c     +                              grdfld( n +  ngx, 4 ) + w * adj( 2 )
c              grdfld( n -  ngx, 4 ) =
c     +                              grdfld( n -  ngx, 4 ) + w * adj( 2 )
c              grdfld( n + ngxy, 4 ) =
c     +                              grdfld( n + ngxy, 4 ) + w * adj( 2 )
c              grdfld( n - ngxy, 4 ) =
c     +                              grdfld( n - ngxy, 4 ) + w * adj( 2 )
c 12 edge centres
c              grdfld( n +  ngx +   1, 4 ) =
c     +                        grdfld( n +  ngx +   1, 4 ) + w * adj( 3 )
c              grdfld( n +  ngx -   1, 4 ) =
c     +                        grdfld( n +  ngx -   1, 4 ) + w * adj( 3 )
c              grdfld( n -  ngx +   1, 4 ) =
c     +                        grdfld( n -  ngx +   1, 4 ) + w * adj( 3 )
c              grdfld( n -  ngx -   1, 4 ) =
c     +                        grdfld( n -  ngx -   1, 4 ) + w * adj( 3 )
c              grdfld( n + ngxy +   1, 4 ) =
c     +                        grdfld( n + ngxy +   1, 4 ) + w * adj( 3 )
c              grdfld( n + ngxy -   1, 4 ) =
c     +                        grdfld( n + ngxy -   1, 4 ) + w * adj( 3 )
c              grdfld( n - ngxy +   1, 4 ) =
c     +                        grdfld( n - ngxy +   1, 4 ) + w * adj( 3 )
c              grdfld( n - ngxy -   1, 4 ) =
c     +                        grdfld( n - ngxy -   1, 4 ) + w * adj( 3 )
c              grdfld( n + ngxy + ngx, 4 ) =
c     +                        grdfld( n + ngxy + ngx, 4 ) + w * adj( 3 )
c              grdfld( n + ngxy - ngx, 4 ) =
c     +                        grdfld( n + ngxy - ngx, 4 ) + w * adj( 3 )
c              grdfld( n - ngxy + ngx, 4 ) =
c     +                        grdfld( n - ngxy + ngx, 4 ) + w * adj( 3 )
c              grdfld( n - ngxy - ngx, 4 ) =
c     +                        grdfld( n - ngxy - ngx, 4 ) + w * adj( 3 )
c 8 corners
c             grdfld( n + ngxy + ngx + 1, 4 ) = grdfld( n + ngxy + ngx + 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n + ngxy + ngx - 1, 4 ) = grdfld( n + ngxy + ngx - 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n + ngxy - ngx + 1, 4 ) = grdfld( n + ngxy - ngx + 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n + ngxy - ngx - 1, 4 ) = grdfld( n + ngxy - ngx - 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n - ngxy + ngx + 1, 4 ) = grdfld( n - ngxy + ngx + 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n - ngxy + ngx - 1, 4 ) = grdfld( n - ngxy + ngx - 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n - ngxy - ngx + 1, 4 ) = grdfld( n - ngxy - ngx + 1, 4 )
c     +                                                    + w * adj( 4 )
c             grdfld( n - ngxy - ngx - 1, 4 ) = grdfld( n - ngxy - ngx - 1, 4 )
c     +                                                    + w * adj( 4 )
           end if
         end do
        end do
      end do
      return
      end
