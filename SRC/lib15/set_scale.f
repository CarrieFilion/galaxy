      subroutine set_scale
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c returns a possible set of conversion factors from my "natural" units to
c   physical units.  Only two can be chosen independently, which this
c   routines assumes will be the length unit and time unit
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      logical gtlogl
c
c local variable
      logical reject
c
c constants in cgs - 4 significant figure precision
      real G, kpc, Msun, Myr
      parameter ( G = 6.672e-8, kpc = 3.086e21,
     +            Msun = 1.989e33, Myr = 3.156e13 )
c
      reject = .true.
      do while ( reject )
        call gtreal( 'Choose unit_L in kpc', unit_L )
        call gtreal( 'Choose unit_T in Myr', unit_T )
c
c calculate unit_V in km/s
        unit_V = unit_L * kpc * 1.e-5 / ( unit_T * Myr )
c
c calculate unit_M in 10^10 solar masses
        unit_M = ( kpc / Myr )**2 * ( kpc * 1.e-10 ) / ( G * Msun )
        unit_M = unit_M * unit_L**3 / unit_T**2
c
        print 200, 'unit_L =', unit_L, ' kpc'
        print 200, 'unit_T =', unit_T, ' Myr'
        print 200, 'unit_V =', unit_V, ' km/s'
        print 200, 'unit_M =', unit_M, ' 10^10 Msun'
        reject = .not. gtlogl( 'Are you happy with these' )
      end do
      scale_set = .true.
      return
 200  format( a, f10.4, a )
      end
