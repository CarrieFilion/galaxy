      subroutine slfset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up tables of self-force and self-PE terms for 2-D and
c   3-D polar grids
c
c See /home/sellwood/p3d/docs/selff.tex for an explanation of these terms
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/supp.f'
c
c externals
      real grnfun, grofu, radfor
c
c local variables
      integer i, ih
      logical done
      real r1, r2
c
      done = .false.
      do ih = 1, ngrid
        call switch( ih )
        if( p2d .or. p3d )then
          if( done )call crash( 'SLFSET', 'Multiple grids' )
c check space
          if( nr( jgrid ) .gt. mrsup )call space(
     +                     mrsup, nr( jgrid ), '/ slffor /', 'SLFSET' )
c work over grid rings
          do i = 1, nr( jgrid )
            r1 = grofu( real( i - 1 ) )
            r2 = grofu( real( i ) )
c self-force terms
            sff( 1, i ) = 2. * pmass * radfor( r1, r1, alpha, 0. )
            sff( 2, i ) = 2. * pmass * ( radfor( r2, r1, alpha, 0. )
     +                                 + radfor( r1, r2, alpha, 0. ) )
            if( p3d )then
              sff( 3, i ) = 4. * pmass * radfor( r1, r1, alpha, dzg )
              sff( 4, i ) = 4. * pmass * ( radfor( r1, r2, alpha, dzg )
     +                                 + radfor( r2, r1, alpha, dzg ) )
            end if
c self-energy terms - not normally used!
            self( 1, i ) = pmass * grnfun( r1, r1, 0., 0. )
            self( 2, i ) = pmass * grnfun( r1, r1, alpha, 0. )
            self( 3, i ) = pmass * grnfun( r2, r1, 0., 0. )
            self( 4, i ) = pmass * grnfun( r2, r1, alpha, 0. )
            if( p3d )then
              self( 7, i ) = pmass * grnfun( r1, r1, 0., dzg )
              self( 8, i ) = pmass * grnfun( r1, r1, alpha, dzg )
              self( 9, i ) = pmass * grnfun( r2, r1, 0., dzg )
              self( 10, i ) = pmass * grnfun( r2, r1, alpha, dzg )
            end if
            if( i .gt. 1 )then
              r2 = grofu( real( i - 2 ) )
              self( 5, i ) = pmass * grnfun( r2, r1, 0., 0. )
              self( 6, i ) = pmass * grnfun( r2, r1, alpha, 0. )
              if( p3d )then
                self( 11, i ) = pmass * grnfun( r2, r1, 0., dzg )
                self( 12, i ) = pmass * grnfun( r2, r1, alpha, dzg )
              end if
            end if
          end do
          done = .true.
        end if
      end do
      if( hybrid )call switch( 0 )
      return
      end
