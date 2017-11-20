      subroutine thracc( jst )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to return an estimate of the central attraction for off-grid
c   particles.  The estimate assumes that the mass distribution still
c   has the anlaytic form of the initial model
c
c calling argument
      integer jst
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c external
      real tsfac
      real*8 frtot, phitot
c
c local variables
      integer i, is
      real fr, r2, tfac, totmas, x, y, z, xp( 3 )
      real*8 rs
      equivalence ( xp( 1 ), x ), ( xp( 2 ), y ), ( xp( 3 ), z )
c
      if( ilist .ne. nlists )call crash ( 'THRACC',
     +                                        'Not off-grid particles' )
      call switch( 0 )
      tfac = tsfac( nzones )
      do is = 1, jst
c get radius
        r2 = 0
        do i = 1, ndimen
          xp( i ) = oldc( i, is ) - xcen( i, jgrid )
          r2 = r2 + xp( i ) * xp( i )
        end do
        rr( is ) = sqrt( r2 )
        if( offthf )then
c compute theoretical radial acceleration in grid units
          rs = rr( is ) / lscale
          fr = frtot( rs ) * ts / gvfac
          gpot( is ) = gpot( is ) + tfac * phitot( rs ) / gvfac**2
        else if( offmnp )then
c use monopole term from on-grid mass only
          call crash( 'THRACC', 'Option not programmed' )
          fr = -totmas / rr( is )**2
          gpot( is ) = gpot( is ) - tfac * totmas / rr( is )
        end if
c resolve into components
        fr = fr / rr( is )
        do i = 1, ndimen
          acc( i, is ) = acc( i, is ) + fr * xp( i )
        end do
c add external perturbation
        if( pertbn )call pertrb( is )
c flag acceleration as complete
        nskip( is ) = .false.
      end do
      return
      end
