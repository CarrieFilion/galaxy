      subroutine scaled
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to multiply the grid mass array or the SFP expansion coeffs
c   by the appropriate consts (pmass) in order to yield values for the
c   acceleration components and potential in grid units.  It is more
c   efficient to do this in one sweep than to normalize the mass of each
c   particle as it is being processed.
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
      include 'inc/model.f'
c
c local variables
      integer i, iuse, j, k, m, n
      real ak, convfac, fac
      save iuse
      data iuse / 0 /
c
      if( ( mzone .le. 0 ) .or. ( mzone .gt. nzones ) )then
        print *, 'mzone =', mzone
        call crash( 'SCALED', 'Impossible mzone' )
      end if
      do jgrid = 1, ngrid
        call switch( jgrid )
c SFP methods
        if( sf2d .or. sf3d )then
c double non-axisymmetric terms because -ve m values are ignored
          if( ( basset .eq. 'bess' ) .or. ( basset .eq. 'lgsp' ) )then
            do i = 1, lastf
              m = msel( i )
              if( m .gt. 0 )then
                n = nsel( i )
                k = n * ngx + maxm + 1
                do j = 1, ngz
                  sfpmss( k - m, 1 ) = 2. * sfpmss( k - m, 1 )
                  sfpmss( k + m, 1 ) = 2. * sfpmss( k + m, 1 )
                  k = k + ngxy
                end do
              end if
            end do
          end if
c rescale to internal units
          convfac = pmass
          if( basset .eq. 'ablj' )convfac = convfac / maxr**2
          do m = -maxm, maxm
            do n = 0, maxn
              fac = convfac
              if( basset .eq. 'bess' )then
                ak = n * deltak
                fac = convfac * sfpwt( n ) * exp( -ak * softl )
              else if( basset .eq. 'lgsp' )then
                fac = convfac * sfpwt( n ) * xcoeff( abs( m ), n, 1 )
              end if
c work over field planes
              j = n * ngx + m + maxm + 1
              do i = 1, ngz
                sfpfld( j, 1 ) = sfpmss( j, 1 ) * fac
                j = j + ngxy
              end do
            end do
          end do
c combine planes if needed
          if( sf3d )call sf3cmb
c PM+SH - scale new coefficients
        else if( s3d )then
          do i = 1, mesh( jgrid )
            s3dfld( i, jgrid ) = pmass * s3dfld( i, jgrid )
          end do
c scale second PM+SH grid for hybrid scheme
          if( hybrid )then
            if(
     +        jgrid .ne. 2 )call crash( 'SCALED', 'hybrid grid mix up' )
            do i = 1, mesh( jgrid )
              s3dfld( i, 1 ) = pmass * s3dfld( i, 1 )
            end do
          end if
c all other grid methods
        else if( lgrd )then
          do i = 1, mesh( jgrid )
            grdmss( i, 1 ) = grdmss( i, 1 ) * pmass
          end do
c Hernquist's field method
        else if( scf )then
          do i = 1, mesh( jgrid )
            sfpfld( i, 1 ) = pmass * sfpmss( i, 1 )
          end do
c direct-N method
        else if( dr3d )then
c set the masses on the first call only
          if( iuse .eq. 0 )then
            icmp = 0
            do i = 1, ncmp
              if( igrd( i ) .eq. jgrid )icmp = i
            end do
            if( icmp .le. 0 )call crash( 'SCALED', 'Pop undefined' )
            if( ndrct .lt. nsp( icmp ) )call crash( 'SCALED',
     +                                          'drpt array too small' )
            do i = 1, nsp( icmp )
              if( uqmass )then
                n = ( i - 1 + kdrct ) * nwpp
                drpt( 1, i ) = pmass * ptcls( n + ncoor + 1 )
              else
                drpt( 1, i ) = pmass
              end if
            end do
            iuse = 1
          end if
        else if( .not. ( noslfg .or. bht ) )then
          call crash( 'SCALED', 'Unrecognized method' )
        end if
      end do
c update accelerations array for field methods
      call switch( 0 )
      return
      end
