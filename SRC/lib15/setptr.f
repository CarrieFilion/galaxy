      subroutine setptr
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to compute the needed sizes of the principal arrays and set
c   some other parameters
c   code extracted from v13 of setrun
      implicit none
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
      integer ih, j
c
c set some defaults
      if( p2d .or. p3d )softl = sqrt( softl2 )
      if( twod )then
        ncoor = 4
      else
        ncoor = 6
      end if
      nplanes = max( nplanes, 1 )
c one list per zone plus a list for offgrid particles
      nlists = ( nzones + nguard ) * ngrid + 1
      if( c3d )then
c separate list for each plane plus one for particles off
        nlists = nplanes + 1
        if( nguard .gt. 0 )call crash( 'SETPTR',
     +                          'intensive care not programmed in c3d' )
        if( ngrid .gt. 1 )then
          if( dr3d )then
            nlists = nlists + 1
          else
            call crash( 'SETPTR',
     +                          'multiple grids not programmed in c3d' )
          end if
        end if
      end if
      if( nlists .gt. lists )then
        if( master )then
          print *, 'needed =', nlists, ' available =', lists
          write( no, * ) 'needed =', nlists, ' available =', lists
        end if
        call crash( 'SETPTR', 'islist array too small' )
      end if
c number of words per particle - coordinates + linked list + population flag
      nwpp = ncoor + 2
      if( uqmass )nwpp = nwpp + 1
      if( bht )nwpp = nwpp + 1
c determine space in / ptcls / needed for particles
      if( numprocs .eq. 1 )then
        lpf = nwpp * nbod
      else
        lpf = nwpp * ( ( nbod - 1 ) / numprocs + 1 )
      end if
c save lptcls
      lptcls = lpf
      if( master )write( no, '( '' ptcls length needed'', i10 )' )lpf
c flag old zones as un-initialized
      if( nzones .gt. 1 )then
        do ih = 1, ngrid
          call switch( ih )
          do j = 2, nzones
            lstep( 1, j ) = -1
            lstep( 2, j ) = -1
          end do
        end do
        call switch( 0 )
      end if
      return
      end
