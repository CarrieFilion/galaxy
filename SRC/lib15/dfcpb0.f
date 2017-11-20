      real*8 function dfcpb0( E, Lz )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c DF for Kalnajs's composite model B0 for the Maclaurin disk
c
c calling arguments
      real*8 E, Lz
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variables
      real*8 arg, Es
      include 'inc/pi.f'
c
c Agris's units:
c   E = .5(r^2 + v_r^2 + v_\theta^2)
c   Lz = r * v_\theta
c where velocities are in units in which \Omega_c = 1 and the potential
c is measured from zero at the centre
c
c My units:
c   Both the potential and the velocities are in units where
c   \Omega_c = \sqrt{ 3\pi / 4 } and the potential is zero at infinity
c
      if( e .ge. 0.5 * Phi0 )then
        if( Lz .gt. 0.d0 )then
          dfcpb0 = 2. * dfnorm( icmp ) / pi
        else
          dfcpb0 = 0
        end if
      else
        Es = E - Phi0
        arg = .75 * pi - 2. * Es + Lz**2
c check for valid argument
        if( arg .le. 0.d0 )then
          print *, E, Lz, Es, Lz**2, arg
          call crash( 'DFCPB0', 'arg .le. 0' )
        end if
        arg = Lz / sqrt( arg )
        dfcpb0 = 2. * dfnorm( icmp ) * ( asin( arg ) + .5 * pi ) / pi**2
      end if
      return
      end
