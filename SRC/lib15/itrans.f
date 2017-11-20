      real*8 function itrans( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine for the integral transform, using the current basis, of the
c   analytic surface density distribution of the initial axisymmetric disk
c called from MASSDI
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/bdpps.f'
c
      common / transi / ak, nt
      integer nt
      real*8 ak
c
c externals
      real*8 ajfun, bsjfun, gsigmt
c
      if( basset .eq. 'bess' )then
c Hankel transform
        itrans = r * bsjfun( 0, r * ak ) * gsigmt( r )
      else if( basset .eq. 'ablj' )then
c Abel-Jacobi transform
        itrans = r * ajfun( 0, nt, r / maxr ) * gsigmt( r )
      else
        call crash( 'ITRANS', 'Unrecognized basis' )
      end if
      return
      end
