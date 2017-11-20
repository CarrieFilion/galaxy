      subroutine codeset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set the logicals in / code / for the appropriate type of
c   field determination method
c as this routine is called from both GETSET and HEDREC, the logicals
c   may have already been set
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local variables
      integer i
c
c set all logicals
      do i = 0, ncodes
        gtype( i ) = .false.
      end do
c work over types
      do i = 1, ngrid
        ncode = igrid( i )
        if( ( ncode .lt. 0 ) .or. ( ncode .gt. ncodes ) )then
          if( lprint .and. master )then
            if( hybrid )then
              print *, ' hybrid type', i, ' set to', ncode
            else if( twogrd )then
              print *, ' twogrd type', i, ' set to', ncode
            else
              print *, 'code type set to', ncode
            end if
          end if
          call crash( 'CODESET',
     +                       'Unrecognised field determination method' )
        end if
        gtype( ncode ) = .true.
c note grid type
        if( lprint .and. master )then
          if( ncode .eq. 1 )print *, '2-D polar grid selected'
          if( ncode .eq. 2 )print *, '3-D Cartesian grid selected'
          if( ncode .eq. 3 )print *, '2-D SFP code selected'
          if( ncode .eq. 4 )print *, '3-D axisymmetric grid selected'
          if( ncode .eq. 5 )print *, '3-D polar grid selected'
          if( ncode .eq. 6 )print *, 'Spherical grid selected'
          if( ncode .eq. 7 )print *, '3-D SFP code selected'
          if( ncode .eq. 8 )print *, '2-D Cartesian grid selected'
          if( ncode .eq. 9 )print *, '3-D SCF code selected'
          if( ncode .eq. 10 )print *, 'direct-N code selected'
          if( ncode .eq. 11 )print *, 'Barnes-Hut tree code selected'
        end if
      end do
c activate all selected methods and set logicals
      call switch( 0 )
      if( ngrid .gt. 1 )then
        if( lsfp )call crash( 'CODESET',
     +                  'Multiple grids cannot include a field method' )
        if( bht )call crash( 'CODESET',
     +                   'Multiple grids cannot include a tree method' )
      end if
c set dimensionality - if not preset
      if( ndimen .le. 1 )then
        twod = c2d .or. p2d .or. sf2d
        if( twod .and. ( hybrid .or. twogrd ) )call crash( 'CODESET',
     +                     'Multiple grids cannot include a 2D method' )
        threed = .not. twod
        if( twod )then
          ndimen = 2
        else
          ndimen = 3
        end if
      else
        if( ndimen .gt. 3 )call crash( 'CODESET',
     +                                          'ndimen > 3 at start!' )
        twod = ndimen .eq. 2
        threed = ndimen .eq. 3
      end if
      ncoor = 2 * ndimen
      return
      end
