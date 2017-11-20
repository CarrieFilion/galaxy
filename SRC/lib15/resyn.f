      real function resyn( t, u )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to reconstruct the disk surface density from its discrete Fourier
c    transform.  It also adds a projected bulge component if bulge = .T.
      use aarrays
      implicit none
c
c calling arguments
      real t, u
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
c externals
      real grofu, sbulge, uofgr
c
c local variables
      integer i, iu, j, m, mm1, n
      real r, th, thfc
c
      if( .not. ( p2d .or. p3d ) )call crash( 'RESYN',
     +                                            'Inappropriate grid' )
      r = grofu( u )
      iu = nint( uofgr( r ) )
      n = ma * iu + 1
      resyn = 0
      if( wres( n ) .ne. 0. )then
c non-axisymmetric part
        m = n + mr * ma
        thfc = alpha * ( t - 1. )
        mm1 = max( mml / nsect, 1 )
        do i = mm1, mmu / nsect
          th = thfc * real( i * nsect )
          j = n + i
          resyn = resyn + wres( j ) * cos( th )
          j = m + i
          resyn = resyn + wres( j ) * sin( th )
        end do
c axisymmetric part
        if( mml .eq. 0 )then
          resyn = wres( n ) * ( 1. + resyn ) / real( na )
          if( bulge )then
            r = grofu( u - .1 ) / lscale
            resyn = resyn + sbulge( r )
          end if
c normalise by undisturbed mass density
c          r = exp( alpha * ( u - 1. ) ) / lscale
c          den = gsigma( r ) * real( na )
c          resyn = resyn * wres( n ) / den
        else
c          resyn = resyn / real( na )
c compute absolute overdensity if requested
          if( .not. bulge )resyn = wres( n ) * resyn ! / real( na )
        end if
      end if
      return
      end
