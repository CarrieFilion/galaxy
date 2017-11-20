      real*8 function velint( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Performs a double integration over all allowed velocities, weighting the
c   DF by preset powers of velocity.  Works for 2-integral DFs in either
c   disk or spherical geometry
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      include 'inc/moment.f'
c
c externals
      real*8 gsigma, phitot, vcirc, velif1
c
c local arrays
      real*8, allocatable :: absc( : ), wt( : )
c
c local variables
      integer i, npts
      real*8 umin, umax
c
      if( ( icmp .le. 0 ) .or. ( icmp .gt. ncmp ) )call crash( 'VELINT',
     +                                        'Nonsense value of icmp' )
c integral is zero outside active disc or if u enters to an odd power
c   or if v enters to an odd power for a non-rotating model
      if( ( r .ge. rmax ) .or. ( mod( iu, 2 ) .ne. 0 ) .or.
     +    ( ( .not. disc( icmp ) ) .and. ( mod( iv, 2 ) .ne. 0 ) ) )then
        velint = 0
c cold disks
      else if( kcold )then
        if( iu .gt. 0 )then
          velint = 0
        else
          velint = gsigma( r )
          if( iv .gt. 0 )velint = velint * vcirc( r )**iv
        end if
      else
c non-adaptive Gauss-Legendre quadrature over allowed half-range of u
        rad = r
        pot = phitot( r )
        umin = 0
        umax = sqrt( 2. * ( phimax - pot ) )
        npts = 16
        allocate ( absc( npts ) )
        allocate ( wt( npts ) )
        call GLqtab( umin, umax, npts, wt, absc )
        velint = 0
        do i = 1, npts
          if( wt( i ) .gt. 0.d0 )velint = velint +
     +                                     wt( i ) * velif1( absc( i ) )
        end do
c double half-range result for disks
        if( disc( icmp ) )velint = 2. * velint
      end if
      return
      end
