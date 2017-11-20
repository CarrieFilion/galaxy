      subroutine rngfrc( rs, rf, zf, ar, az, pot )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the radial and vertical components of the acceleration and the
c   potential at the field point (rf,zf) due to a unit mass ring at the
c   source point (rs,0).  The elliptic integral expressions are similar
c   to those given in Binney & Tremaine p72
c
c calling arguments
      real*8 ar, az, pot, rs, rf, zf
c
c externals - complete elliptic integrals of the first and second kind
c  these use the Numerical Recipes routine which is faster than NAG
      real*8 celli1, celli2
c
      real*8 c1, c2, der, hm, k, k2, zt
      include 'inc/pi.f'
c
c special values
      if( rs .eq. 0. )then
        k = rf**2 + zf**2
        if( k .eq. 0. )then
          ar = 0
          az = 0
          pot = -1.
        else
          k = 1. / sqrt( k )
          pot = -k
          k = k**3
          ar = -rf * k
          az = -zf * k
        end if
      else
c avoid divergent result if the field point coincides with the source
        zt = zf
        if( ( zf .eq. 0. ) .and. ( rf .eq. rs ) )zt = .01 * rs
c radial force
        if( rf .gt. 0. )then
          k2 = 4. * rf * rs / ( ( rs + rf )**2 + zt**2 )
          k = sqrt( k2 )
          c1 = k * celli1( k )
          c2 = celli2( k ) / ( 1. - k2 )
          der = .125 * k**3 *
     +                       ( ( rs + zt**2 / rs ) / rf - rf / rs ) / rf
          hm = sqrt( rs * rf )
          ar = ( c2 * der - .5 * c1 / rf ) / ( pi * hm )
          pot = -c1 / ( pi * hm )
c z force - zero value in symmetry plane
          if( zf .eq. 0. )then
            az = 0
          else
            az = -c2 * k**3 * zf / ( 4. * pi * hm**3 )
          end if
        else
c values on symmetry axis
          ar = 0
          k = 1. / sqrt( rs**2 + zf**2 )
          az = -zf * k**3
          pot = -k
        end if
      end if
      return
      end
