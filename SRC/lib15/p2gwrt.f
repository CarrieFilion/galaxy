      subroutine p2gwrt( azimth, lwork, w )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c routine to write out transformed Green function for 2-D poalr grid
c   The Fourier coefficients were calculated in p2grnm and stored in
c   / ptcls /  This routine picks them out and writes them in the order
c   which will be most useful for the convolution
c
c calling arguments
      integer lwork
      logical azimth
      real w( lwork )
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
      integer ig, ip, ir, j, k, m, n, nm, nrg
c
      nrg = nr( jgrid )
      nm = ng
      if( azimth )nm = min( mmax - 2, ng - 1 )
      k = nrg * nm
      if( k .gt. lwork )call space( lwork, k, 'w', 'P2GWRT' )
c scan over sectoral harmonics
      do 1 ig = 1, ng
      if( ( ( ig .eq. 1 ) .or. ( ig .eq. mmax ) ) .and. azimth )go to 1
c build matrix nr * nr
        n = 1 - k
        do ir = 1, nrg
          n = n + k
          call blkcpy( ptcls( n ), w( 1 ), k )
          j = ir - nrg
          m = ig - nm
          if( azimth )m = m - 1
          do ip = 1, nrg
            j = j + nrg
            m = m + nm
            grdmss( j, 1 ) = w( m )
          end do
        end do
c write matrix to file
        write( ngt )( grdmss( m, 1 ), m = 1, j )
    1 continue
      return
      end
