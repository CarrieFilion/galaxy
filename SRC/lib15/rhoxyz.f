      real*8 function rhoxyz( x, y, z )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Evaluates the volume density at an arbitrary point.
c   The position coordinates and returned value are in natural units.
c
c calling arguments
      real*8 x, y, z
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 rhohal, rhokks
c
c local variables
      integer i
      real*8 r
c
c Kuz'min-Kutuzov oblate isochrone
      if( ctype( icmp ) .eq. 'KKSP' )then
        r = sqrt( x**2 + y**2 )
        rhoxyz = rhokks( r, z )
      else
c use generic function
        r = sqrt( x**2 + y**2 + z**2 )
        rhoxyz = 0
        i = icmp
        do icmp = 1, ncmp
          if( .not. disc( icmp ) )rhoxyz = rhoxyz + rhohal( r )
        end do
        icmp = i
      end if
      return
      end
