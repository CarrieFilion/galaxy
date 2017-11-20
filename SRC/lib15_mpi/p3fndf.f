      subroutine p3fndf( jzone, kzone )
c  Copyright (C) 2015, Jerry Sellwood
c
c Determines the forces and potentials arising from a single zone of masses
c   by Fourier transformation on the 3-D polar grid
c The vertical convolution requires the mesh to be doubled in order to remove
c   the effects of the periodic images.  The extra space is allocated between
c   the potential array and the mass arrays, so the mass array has to be
c   copied there before calling P3MANL.
c The routine calls P3MANL to transform and double the size of the mass
c   array.  It then copies the array to each of the acceleration component
c   arrays for the radial part of the convolution (P3CONV) and resynthesis
c   by P3FSYN.
c As the arrays for the radial, azimuthal, vertical acceleration components
c   and also the potentials are contiguous, the extra space needed during
c   the solution for each can be borrowed from the next array.
      use aarrays
      implicit none
c
c calling arguments
      integer jzone, kzone
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
      integer i, ips, ispac, itype, j, jrf, jrs, krf, krs, ntypes
      logical skip
c
c determine number of active sectoral harmonics
      j = 0
      do i = 1, ng
        if( .not. lg( i ) )j = j + 1
      end do
      ntypes = 3
      if( potl )ntypes = 4
      skip = j .eq. 0
      if( .not. skip )then
c determine radial range which contains all sources
        jrs = nr( jgrid )
        krs = 0
        do i = 1, nzones
          jrs = min( jrs, jrad( i ) )
          krs = max( krs, krad( i ) )
          if( i .gt. 1 )krad( i ) = max( krad( i ), krs )
        end do
        if( ( jrs .le. krs ) .and. ( krs .gt. 0 ) )then
c copy mass array to work space
          ispac = ( krs - jrs + 1 ) * na * ngz
          ips = ngz * na * ( jrs - 1 )
          call blkcpy( grdmss( ips + 1, 1 ),
     +                 grdmss( 1, 0 ), ispac )
c form double Fourier transform of mass array
          if( parallel )then
            call p3manlp( jrs, krs )
          else
            call p3manl(jrs, krs )
          end if
c assume forces are required over the whole grid unless told otherwise
          jrf = 1
          if( ( jzone .gt. 0 ) .and.
     +        ( jrad( jzone ) .lt. nr( jgrid ) ) )jrf = jrad( jzone )
          krf = nr( jgrid )
          if( ( kzone .le. nzones ) .and.
     +        ( .not. wholeg ) )then
            krf = krad( kzone )
            do while ( krf .lt. jrf )
              kzone = kzone + 1
              if(
     +       kzone .gt. nzones )call crash( 'P3FNDF', 'kzone > nzones' )
              krf = krad( kzone )
            end do
          end if
c work over types - vertical forces only if fixrad is .true.
          itype = 0
          if( fixrad )itype = 2
          do while ( itype .lt. ntypes )
            itype = itype + 1
c no azimuthal forces when m = 0 only is selected
            if( ( itype .eq. 2 ) .and. ( ng .eq. 1 ) )then
              do i = 1, mesh( jgrid )
                grdfld( i, 2 ) = 0
              end do
            else
c convolve with Green function
              if( parallel )then
                call p3convp( itype, jrs, krs, jrf, krf )
              else
                call p3conv( itype, jrs, krs, jrf, krf )
              end if
c re-synthesize
              call p3fsyn( itype, jrf, krf )
c this call not implemented as the parallel version is slower
c              call p3fsynp( itype, jrf, krf, lmesh, msh )
            end if
          end do
        else if( krs .gt. nr( jgrid ) )then
          call crash( 'P3FNDF', 'Mass outside the grid?' )
        else
          if( master )write( no, * )' No mass assigned to grid', istep
          skip = .true.
c          call
c     +       crash( 'P3FNDF', 'Invalid radial range for source masses' )
        end if
      end if
c mass-free or force-free grid
      if( skip )then
        do itype = 1, ntypes
          do i = 1, mesh( jgrid )
            grdfld( i, itype ) = 0
          end do
        end do
      end if
      return
      end
