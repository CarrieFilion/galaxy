      real function tsfac( iz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns time step factor, relative to the standard step size, for
c   the given zone
c
c calling argument
      integer iz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local variable
      integer i
c
c particle in guard zone
      if( iz .le. 0 )then
        i = 1 - iz
        if( i .gt. nguard )call crash( 'TSFAC',
     +                                  'Impossible guard zone number' )
        tsfac = 1. / real( nsub( i ) )
      else
        if( iz .gt. nlists )call crash( 'TSFAC',
     +                                        'Impossible zone number' )
c no time step zones allowed in c3d
        if( c3d )then
          tsfac = 1
        else
          i = min( iz, nzones )
          tsfac = nstep( i )
        end if
      end if
      return
      end
