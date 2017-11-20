      real*8 function phikks( r, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns gravitational potential of a flattened Kuz'min-Kutuzov oblate
c   spheroid.   Calling arguments are cylindrical radius and z.  The value
c   of a is hrad and the axis ratio c/a is haloc; both are stored in / model /
c   The notation is that used by Dejonghe & de Zeeuw in Ap J v333 p90 (1988).
c
c calling arguments
      real*8 r, z
c
c local variables
      real*8 lambda, nu
c
      call kkcoor( r, z, lambda, nu )
c
      phikks = -1. / ( sqrt( lambda ) + sqrt( nu ) )
      return
      end
