      subroutine findf( lpot )
c  Copyright (C) 2015, Jerry Sellwood
c
c Driver routine to compute the forces and potential for a generic method.
c   The calling argument flags whether the potential is, or is not, needed.
c This routine does almost nothing when the field is determined by an SFP or
c   PM+SH method
      implicit none
c
c calling argument
      logical lpot
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
      integer i, jzone, kzone
      logical full
c
      if( ( mzone .lt. 1 ) .or. ( mzone .gt. nzones ) )then
        print *, 'mzone =', mzone
        call crash( 'FINDF', 'invalid mzone' )
      end if
c save calling argument
      potl = lpot
c
      do i = 1, ngrid
        call switch( i )
        if( ncode .gt. 11 )call crash( 'FINDF', 'Unrecognised method' )
c skip if nothing to do
        if( .not. ( s3d .or. lsfp ) )then
c polar grids
          if( p2d .or. p3d )then
c set flags for full or part grid calculation and report if requested
            full = potl .or. ( mzone .eq. nzones ) .or. wholeg
            if( full )then
              jzone = 0
              kzone = nzones + 1
              if( lprint .and. master )write( no,
     +                        '( '' Full grid force determination'' )' )
            else
              jzone = 1
              kzone = mzone
              if( lprint .and. master )write( no,
     +            '( '' Force determination: partial grid from'', i3,'//
     +                     ' '' to'', i3)' )jrad( jzone ), krad( kzone )
            end if
c compute fields
            if( p2d )then
              call p2fndf( jzone, kzone )
            else
              call p3fndf( jzone, kzone )
            end if
c convert to Cartesian components
            call polcat
c 3D Cartesian grid
          else if( c3d )then
            call c3fndf
c axisymmetric polar grid
          else if( p3a )then
            call pafndf
c 2D Cartesian grid
          else if( c2d )then
            call c2fndf
c direct-N
          else if( dr3d )then
            call drctN2
c Barnes-Hut tree method
          else if( bht )then
            call bhtree
          end if
        end if
      end do
c report if requested
      if( lprint .and. master .and. .not. ( s3d .or. lsfp ) )then
        if( .not. potl )then
          write( no, '( '' Found forces at step'', i7 )' )istep
        else
          write( no,
     +   '( '' Found forces and potentials at step'', i7 )' )istep
        end if
      end if
c set flag that gravitational fields are now current
      isfld = istep
      call switch( 0 )
      return
      end
