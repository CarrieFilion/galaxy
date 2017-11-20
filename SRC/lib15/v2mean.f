      real*8 function v2mean( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the mean-square azimuthal velocity in a disk or sphere at the
c   requested radius
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 phitot, sigmau, vcirc, vmean
c
c local variables
      integer i, n
      real*8 a, dfm, r1, r2
c
      dfm = dfcns( 1, icmp )
c cold disk
      if( kcold )then
        v2mean = vcirc( r )**2
c Kalnajs function
      else if( cdft( icmp ) .eq. 'KALN' )then
c Kalnajs's equation (23) can be rewritten as $$ \eqalign{
c r\Sigma(r)\sigma_u^2(r) & = -{1 \over r^m} \int_0^r r^{\prime m}
c \Sigma(r^\prime) \left[ \phi(r^\prime) + V_c^2(r^\prime) \right] dr^\prime \cr
c & = -{1 \over r^m} \int_0^r f(r^\prime) dr^\prime. \cr}
c $$ Thus $$ \eqalign{
c {d \over dr}(r\Sigma\sigma_u^2) & = {m \over r^{m+1}} \int_0^r f(r^\prime)
c dr^\prime - {f(r) \over r^m} \cr
c & = -m \Sigma\sigma_u^2 - \Sigma \left[ \phi + V_c^2 \right] \cr}
c $$ and $$
c {1 \over \Sigma} {d \over dr}(r\Sigma\sigma_u^2) = -m\sigma_u^2-\phi-V_c^2.
c $$  Now, the Jeans equation gives us $$
c {1 \over \Sigma} {d \over dr}(r\Sigma\sigma_u^2) = \langle v^2 \rangle - V_c^2
c $$ so therefore $$
c \langle v^2 \rangle = -m\sigma_u^2 - \phi.
c $$
        v2mean = -dfm * sigmau( r )**2 - phitot( r )
c Miyamoto functions for KT disk - formula (33)
      else if( cdft( icmp ) .eq. 'MIYA' )then
        v2mean = 1. / ( dble( 2 * mmiy + 4 ) * sqrt( 1. + r * r ) ) +
     +  dble( mmiy ) * r**2 / ( dble( mmiy + 2 ) * ( 1. + r * r )**1.5 )
c Omega model - vely distribution is isotropic in plane
      else if( cdft( icmp ) .eq. 'OMEG' )then
        v2mean = vmean( r )**2 + sigmau( r )**2
c Evans's DFs for power law discs
      else if( cdft( icmp ) .eq. 'EVAN' )then
        v2mean = ( 1. + dfm ) /
     +                        ( ( 2. * evbeta + dfm + 1. ) * r**evbeta )
c Evans/Collett DFs for Rybicki discs - MN v264 p353
      else if( cdft( icmp ) .eq. 'EVCO' )then
        n = nint( dfm )
        r2 = r * r
        r1 = sqrt( 1. + r2 )
c Appendix C - the result from this loop does not look quite right!
        if( n .ge. 0 )then
          v2mean = 1
          a = 2**n * r1**3 / ( 1. + r1 )**n
          do i = 0, n
            if( i .gt. 0 )a = a * r2 * dble( n + 1 - i )
     +                                 / ( ( 1. + r1 ) * dble( 2 * i ) )
            v2mean = v2mean + a * dble( 2 * i + 1 ) / dble( n + i + 1 )
          end do
          v2mean = v2mean / ( ( 1. + r2 ) * ( 1. + r1 ) )
        else
c equation (2.8)
          v2mean = r1 * log( ( 1. + r1 ) / r1 )
        end if
      else
        call crash( 'V2MEAN', 'Unknown DF type' )
      end if
      return
      end
