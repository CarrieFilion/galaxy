      subroutine p3gwrt( itype )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to write out transformed Green function for 3-D polar grid
c   The Fourier coefficients were calculated in P3GRNM and stored in
c   / ptcls /  This routine picks them out and writes them in the order
c   which will be needed for P3CONV
c the ouput file has already been opened and the header written in P3GRNM
c
c calling argument
      integer itype
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c local variables
      integer im, im1, im2, ir, is, iz, iz1, iz2, j, k, m, n, nm, nz
      logical azimth, vert
c
      azimth = itype .eq. 2
      vert = itype .eq. 3
c
      im1 = 1
      im2 = ng
      if( azimth )then
        im1 = 2
        im2 = min( ng, mmax - 1 )
      end if
      nm = im2 - im1 + 1
c
      iz1 = 1
      iz2 = ngz + 1
      if( vert )then
        iz1 = 2
        iz2 = ngz
      end if
      nz = iz2 - iz1 + 1
      n = nm * nz
c work over z
      do iz = iz1, iz2
c work over m
        do im = 1, nm
c build matrix nr * nr
          j = 0
          k = nm * ( iz - iz1 ) + im
          do is = 1, nr( jgrid )
            do ir = 1, nr( jgrid )
              j = j + 1
              grdmss( j, 1 ) = ptcls( k )
              k = k + n
            end do
          end do
c write matrix to file
          write( ngt )( grdmss( m, 1 ), m = 1, j )
        end do
      end do
      return
      end
