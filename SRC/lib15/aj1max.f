      real*8 function aj1max( aj2 )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes the maximum possible radial action for the given angular
c   momentum (aj2)
c
c calling argument
      real*8 aj2
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 actj1, Emax
c
c local variables
      real*8 E
c
c check whether this angular momentum is allowed
      if( abs( aj2 ) .gt. Lztmx( icmp ) )then
        aj1max = 0
      else
c compute radial action along energy cut off
        E = Emax( aj2 )
        aj1max = actj1( E, aj2 )
      end if
      return
      end
