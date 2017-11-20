      subroutine getacc( jst )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c Returns total accelerations acting on each particle in the current group.
c Calls SFPF, GRIDF or S3DFOR to determine the self-consistent part of the field
c   and any supplementary or external components are added in this routine.
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c externals
      real axifrc, frcorr, halfrc, tsfac
c
c local variables
      integer i, ig, is, n
      logical lshort
      real fr, rad, t, x, y, z, xp( 3 )
      equivalence ( xp( 1 ), x ), ( xp( 2 ), y ), ( xp( 3 ), z )
c
      if( istep .ne. isfld )then
        print *, istep, isfld
        call crash( 'GETACC', 'Accelerations at useless time' )
      end if
c clear the acceleration and potential arrays
      lshort = .false.
      do is = 1, jst
        do i = 1, 3
          acc( i, is ) = 0
        end do
        gpot( is ) = 0
c raise flag if any of these particles is in the guard zone
        lshort = lshort .or. ( iz( is ) .le. 0 )
      end do
c work over active methods
      do ig = 1, ngrid
        call switch( ig )
        if( ncode .gt. 11 )call crash( 'GETACC', 'Unrecognised method' )
c flag all particles as active - forces come from rigid components
        if( noslfg )then
          do is = 1, jst
            nskip( is ) = .true.
          end do
c tree method
        else if( bht )then
          call bhaccl( jst )
c direct-N method
        else if( dr3d )then
          if( stphev )then
c copy accelerations and potential from precalculated array
            do is = 1, jst
              n = 1 + loc( is ) / nwpp - kdrct
              do i = 1, 3
                acc( i, is ) = acc( i, is ) + drpt( i + 7, n )
              end do
              gpot( is ) = gpot( is ) + drpt( 11, n )
              nskip( is ) = .true.
            end do
          else
c add forces from heavies
            call hevacc( jst )
          end if
c skip grid forces on heavy particles - two pops interact directly
        else if( .not. stphev )then
c SFP methods
          if( sf2d .or. sf3d )then
            call sfpacc( jst )
c SCF method
          else if( scf )then
            call scfacc( jst )
c PM+SH method
          else if( s3d)then
c$$$            if( izone .gt. 0 )then
c$$$              call s3dacc( jst )
c$$$            else
              call s3dacc0( jst )
c$$$            end if
c grid methods
          else
            call grdacc( jst )
          end if
        end if
      end do
c restore currently active grid
      call switch( jlist )
c
c sum forces between on- and off-grid particles if requested
      if( offfor )then
        call offacc( jst )
c various force rules for off-grid particles - no forces if no rule (default)
      else if( ilist .eq. nlists )then
        if( offthf .or. offmnp )call thracc( jst )
      end if
c
      do is = 1, jst
        if( nskip( is ) )then
          rr( is ) = 0
          do i = 1, ndimen
            xp( i ) = oldc( i, is ) - xcen( i, jgrid )
            rr( is ) = rr( is ) + xp( i )**2
          end do
          rr( is ) = sqrt( rr( is ) )
c fixed forces
          fr = 0
c supplementary forces (if any) - assumed spherical
          if( fixrad )then
            fr = axifrc( rr( is ) )
          else
            if( rigidh )fr = halfrc( rr( is ) )
            if( suppl )fr = fr + frcorr( rr( is ) )
          end if
          if( fr .ne. 0. )then
            fr = fr / rr( is )
            do i = 1, ndimen
              acc( i, is ) = acc( i, is ) + fr * xp( i )
            end do
          else if( fixrad )then
c cylindrically symmetric forces
            rad = sqrt( x * x + y * y )
            fr = axifrc( rad ) / rad
            do i = 1, 2
              acc( i, is ) = acc( i, is ) + fr * xp( i )
            end do
          end if
c add external perturbation
          if( pertbn )call pertrb( is )
c end skip of particles off grid
        end if
      end do
c rescale accelerations for zone
      if( nzones .gt. 1 )then
        do is = 1, jst
          if( iz( is ) .gt. 1 )then
            t = tsfac( iz( is ) )**2
            do i = 1, ndimen
              acc( i, is ) = t * acc( i, is )
            end do
          end if
        end do
      end if
c particles on extra-short time steps
      if( lshort )then
        do is = 1, jst
          if( iz( is ) .le. 0 )then
            t = tsfac( iz( is ) )**2
            do i = 1, ndimen
              acc( i, is ) = t * acc( i, is )
            end do
          end if
        end do
      end if
      return
      end
