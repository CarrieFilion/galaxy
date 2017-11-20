      subroutine densty
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to write information about the current density distribution into
c   the .res file.  The data written are the expansion coefficients for SFP
c   methods, or the density array or its Fourier transform for grids.
c For models with more than a single time step zone, or for 3-D models,
c   the density has to be summed over the grids, or planes
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      common / rjmass / factor, scfact
      real factor, scfact
c      include 'inc/rjmass.f'
c
c local allocatable array
      real, allocatable :: w(:)
c
c external
      real grofu
c
c local variables
      character*4 bstr( 6 )
      logical lfixrad
      integer i, ih, iz, j, k, l, m, n
      real ahol, area, c, convfac, r1, r2, sigmaf, xx
      include 'inc/pi.f'
c
      data bstr / 'PROJ', 'DANL', 'DNST', 'SFPC', 'SCFC', 'S3DC' /
c
      l = 0
      do ih = 1, ngrid
        l = max( l, mesh( ih ) )
      end do
      allocate ( w( l ) )
      do ih = 1, ngrid
        call switch( ih )
c SFP methods
        if( sf2d .or. sf3d )then
          convfac = 1. / ( lscale**3 * ts**2 )
c sum over planes
          do k = 1, 2 * lastf
            w( k ) = 0
          end do
          do j = 1, ngz
            k = 0
            do i = 1, lastf
              m = msel( i )
              n = nsel( i )
              l = ( j - 1 ) * ngxy + n * ngx + maxm + 1
c m=0 terms are done twice - necessary for consistency with previous versions
              w( k + 1 ) = w( k + 1 ) + convfac * sfpmss( l + m, 1 )
              w( k + 2 ) = w( k + 2 ) + convfac * sfpmss( l - m, 1 )
              k = k + 2
            end do
          end do
          if( master )then
            read( bstr( 4 ), '( a4 )' )ahol
            write( nphys )irun, ahol, istep, basis, lastf
            write( nphys )maxr, ( w( i ), i = 1, 2 * lastf ), ncontrib
          end if
c SCF method
        else if( scf )then
          convfac = 1. / ( lscale**3 * ts**2 )
          do i = 1, mesh( jgrid )
            w( i ) = convfac * sfpmss( i, 1 )
          end do
          if( master )then
            read( bstr( 5 ), '( a4 )' )ahol
            write( nphys )irun, ahol, istep, basis, mesh( jgrid )
            write( nphys )( w( i ), i = 1, mesh( jgrid ) ), ncontrib
          end if
c S3D method
        else if( s3d .and. s3dc )then
          convfac = 1. / ( lscale**3 * ts**2 )
          do i = 1, mesh( jgrid )
            w( i ) = convfac * s3dfld( i, 1 )
          end do
          if( master )then
            read( bstr( 6 ), '( a4 )' )ahol
            write( nphys )irun, ahol, istep, nr( jgrid ), s3lmax
            write( nphys )( w( i ), i = 1, mesh( jgrid ) )
          end if
        end if
        if( dnst .or. danl )then
c set constant for surface density conversion
          sigmaf = 1. / ( lscale * ts * ts )
c convert mass distribution to surface density
          if( p2d .or. p3d )then
            l = 0
            r1 = grofu( 0. )
            do i = 1, nr( jgrid )
              r2 = i
              r2 = grofu( r2 - .5 )
c .5 * alpha = pi / na
              area = .5 * alpha * ( r2 * r2 - r1 * r1 )
              c = sigmaf / area
              r1 = r2
              do j = 1, na
                l = l + 1
                if( p2d )then
                  xx = grdmss( l, 1 )
                else
c project onto a single plane
                  n = j + ( i - 1 ) * na * ngz
                  xx = 0.
                  do iz = 1, ngz
                    xx = xx + grdmss( n, 1 )
                    n = n + na
                  end do
                end if
c borrow space from radial forces
                grdfld( l, 1 ) = c * xx
              end do
            end do
c write out raw density array
            if( dnst )then
              if( master )then
                read( bstr( 3 ), '( a4 )' )ahol
                write( nphys )irun, ahol, istep, nr( jgrid ), na
                k = nr( jgrid ) * na
                write( nphys )( grdfld( i, 1 ), i = 1, k )
              end if
            end if
c project mass distribution onto the three grid planes
          else if( c3d )then
            if( dnst )then
c borrow space from potential array
              n = ngy * ngx + ngz * ngx + ngz * ngy
              do i = 1, n
                grdfld( i, 4 ) = 0
              end do
c x-y projection
              k = 0
              do j = 1, ngz
                do i = 1, ngxy
                  k = k + 1
                  grdfld( i, 4 ) = grdfld( i, 4 ) + grdmss( k, 1 )
                end do
              end do
c x-z projection
              k = 0
              l = ngy * ngx
              do j = 1, ngz * ngy
                do i = 1, ngx
                  l = l + 1
                  k = k + 1
                  grdfld( l, 4 ) = grdfld( l, 4 ) + grdmss( k, 1 )
                end do
                if( mod( j, ngy ) .ne. 0 )l = l - ngx
              end do
c y-z projection
              k = 0
              l = ngy * ngx + ngz * ngx
              do j = 1, ngz * ngy
                l = l + 1
                do i = 1, ngx
                  k = k + 1
                  grdfld( l, 4 ) = grdfld( l, 4 ) + grdmss( k, 1 )
                end do
              end do
c normalise and save
              sigmaf = sigmaf * factor / sqrt( 6. )
              do i = 1, n
                grdfld( i, 4 ) = sigmaf * grdfld( i, 4 )
              end do
              i = 3
              if( master )then
                read( bstr( 1 ), '( a4 )' )ahol
                write( nphys )irun, ahol, istep, i, n
                write( nphys )( grdfld( i, 4 ), i = 1, n )
              end if
            end if
c 3-D axisymmetric code
          else if( p3a )then
c sum over planes in z - borrow space from radial forces
            l = 0
            r1 = grofu( 0. )
            do i = 1, nr( jgrid )
              r2 = i
              r2 = grofu( r2 - .5 )
              area = pi * ( r2 * r2 - r1 * r1 )
              c = sigmaf / area
              r1 = r2
              xx = 0.
              k = i - nr( jgrid )
              do j = 1, ngz
                k = k + nr( jgrid )
                xx = xx + grdmss( k, 1 )
              end do
              l = l + 1
              grdfld( l, 1 ) = c * xx
            end do
c write out raw density array
            if( dnst )then
              i = 1
              if( master )then
                read( bstr( 3 ), '( a4 )' )ahol
                write( nphys )irun, ahol, istep, nr( jgrid ), i
                write( nphys )( grdfld( i, 1 ), i = 1, nr( jgrid ) )
              end if
            end if
          end if
        end if
        if( danl .and. ( p2d .or. p3d ) )then
c save the mass array and replace it with the scaled density array
          n = na * nr( jgrid )
          call blkcpy( grdmss( 1, 1 ), grdfld( 1, 2 ), n )
          call blkcpy( grdfld( 1, 1 ), grdmss( 1, 1 ), n )
c find Fourier coefficients
c use 2-D routine even for p3d code since mass is projected
          lfixrad = lg( 1 )
          lg( 1 ) = .false.
c direct data ready in temporary space
          call p2manl( 1, 1, nr( jgrid )  )
          lg( 1 ) = lfixrad
c rearrange and normalize
          l = 0
          n = ng * nr( jgrid )
          do i = 1, nr( jgrid )
            k = ( i - 1 ) * na + 1
c m=0 term
            c = grdmss( k, 1 )
            l = l + 1
            n = n + 1
            w( l ) = c
            w( n ) = 0.
c active Fourier components
            do j = 2, na, 2
              m = j / 2 + 1
              if( m .le. ng )then
                l = l + 1
                n = n + 1
                w( l ) = 0.
                w( n ) = 0.
                if( c .gt. 0. )then
                  xx = max( abs( grdmss( k + 1, 1 ) ),
     +                      abs( grdmss( k + 2, 1 ) ) )
                  if( c .gt. 1.e-6 * xx )then
                    w( l ) = grdmss( k + 1, 1 ) / c
                    w( n ) = -grdmss( k + 2, 1 ) / c
                  end if
                end if
              end if
              k = k + 2
            end do
          end do
c write out Fourier coeffs
          if( master )then
            read( bstr( 2 ), '( a4 )' )ahol
            write( nphys )irun, ahol, istep, nr( jgrid ), ng
            write( nphys )( w( i ), i = 1, n )
          end if
c restore mass array
          n = na * nr( jgrid )
          call blkcpy( grdfld( 1, 2 ), grdmss( 1, 1 ), n )
        end if
      end do
      call switch( 0 )
c return workspage
      deallocate ( w )
      return
      end
