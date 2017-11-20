      real*8 function vprob( vmax )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns a random v/sigma selected from a Maxwellian velocity distribution
c   with upper bound vmax
c   uses QUADPACK routines DQAGS and DQNG
c
c calling argument
      real*8 vmax
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c externals
      external vmxwel
      real*8 quad_gnr, quad_nad, ranuni
c
c local array
      real*8 c( 4 )
c
c local variables
      integer ier, ind, ir, jfail
      real*8 epsa, epsr, err, f, fmax, reqf, tol, v1, v2, vs, zero
      parameter ( zero = 0 )
      save vs, fmax
      data vs / 0.d0 /
c
      if( vs .ne. vmax )then
c integrate Maxwellian to set new max value
        epsa = 1.e-6
        epsr = 1.e-6
        fmax = quad_gnr( vmxwel, zero, vmax, epsa, epsr, ier )
        if( ier .ne. 0 )then
          print *, 'ier =', ier, ' from QUAD_GNR'
          call crash( 'VPROB', 'QUADPACK error' )
        end if
        vs = vmax
      end if
c select random fraction
      reqf = fmax * ranuni( 0. )
c find velocity for this fraction
      v1 = 0
      v2 = vmax
      tol = 1.e-6
      ind = 1
      ir = 0
      jfail = 0
      do while ( ind .ne. 0 )
        call fndzro( v1, v2, err, tol, ir, c, ind, jfail )
c use a slightly faster, but non adaptive, integrator
        f = quad_nad( vmxwel, zero, v1, epsa, epsr, ier )
c ier = 1 means requested accuracy was not achieved
        if( ier .gt. 1 )then
          print *, 'ier =', ier, ' from QUAD_NAD'
          call crash( 'VPROB', 'QUADPACK error' )
        end if
        err = f - reqf
      end do
      vprob = v1
      return
      end
