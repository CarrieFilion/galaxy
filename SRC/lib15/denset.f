      subroutine denset( norm )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to set up disc and halo densities - called from dfiter
c
c the disc mass distribution is calculated once and stored separately.
c   The halo mass distribution is set from an integration of the DF over
c   all velocities.  The mass ratio is adjusted here after halo assignment
c   and the DF normalization is changed in order to converge to the desired
c   mass ratio.
      use aarrays
      implicit none
c
c calling argument
      real norm
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/dfitcm.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c local allocatable array
      integer lmeshd
      real, allocatable :: mshd(:)
      save mshd
c
c externals
      external dskrho, sphrho
      real*8 distfn, gmassi, Phispl
c
c local variables
      integer i, icall, j, jcmp, jdisc, kgrid, l, n
      logical first
      real b1, b2, m0, m1, m2, totdisc, tothalo
      real*8 edge, zero
      parameter ( zero = 0. )
      save b1, b2, first, icall, jdisc, m1, m2, totdisc
c
      data first / .true. /
c
      kgrid = jgrid
      if( first )then
        lmeshd = mesh( jgrid )
        allocate ( mshd( lmeshd ) )
        totdisc = 0
c select the disk component
        jdisc = -1
        do j = 1, ncmp
          if( disc( j ) )jdisc = j
        end do
        if( jdisc .lt. 0 )call crash( 'DENSET', 'No disk' )
        print *, 'Mass component', jdisc, ' selected for grid',
     +                                  jgrid, ', disc =', disc( jdisc )
        jcmp = icmp
        icmp = jdisc
        kcold = disc( icmp ) .and. cdft( icmp ) .eq. 'NONE'
        if( s3d .and. ( nr( jgrid ) * s3ntm * 4 .gt. lmeshd )
     +                   )call crash( 'DENSET', 'Temp array too small' )
        rmax = rtrunc( icmp )
        edge = rgrid( jgrid ) / lscale
        if( rmax .gt. edge )call crash( 'DENSET',
     +                                         'disc larger than grid' )
c set mass arrays to zero
        call masscl( .false. )
c assign disc mass to grid in appropriate manner for field method
        j = igrd( jdisc )
        if( hybrid )then
          j = 2
          lhyb2 = .true.
        end if
        call switch( j )
        if( c3d )then
          call massdf( dskrho )
        else if( s3d )then
          call mass3d
c save disc mass array
          do i = 1, mesh( jgrid )
            mshd( i ) = s3dmss( i, 1, jgrid )
          end do
        else
          call massdi
        end if
        call mascmb
        call switch( kgrid )
c get total disc mass
        if( sf3d )then
          totdisc = gmassi( rmax )
        else if( s3d )then
          n = nr( jgrid )
          j = 0
          l = 4 * s3ntm * ( n - 1 ) + 1
          totdisc = s3dfld( 1, jgrid ) + s3dfld( l, jgrid ) * s3rad( n )
        else if( lgrd )then
          totdisc = 0.
          do i = 1, mesh( jgrid )
            mshd( i ) = grdmss( i, 1 )
            totdisc = totdisc + grdmss( i, 1 )
          end do
        else
          call crash( 'DENSET', 'Unrecognized option' )
        end if
c mass array was in program units
        if( .not. sf3d )totdisc = totdisc / ( lscale**3 * ts**2 )
        print * , 'disc mass assigned = ', totdisc
c force active mass ratio to be equal to requested mass ratio
        fmass( jdisc ) = totdisc / cmpmas( jdisc )
        fmass( jcmp ) = fmass( jdisc )
        first = .false.
        icall = 0
        icmp = jcmp
        call switch( kgrid )
c check outer edge
        edge = min( rgrid( jgrid ), zm( jgrid ) ) / lscale
        if( edge .lt. rtrunc( icmp ) )then
          print *, 'Warning, component', icmp, ' extends to',
     +                         rtrunc( icmp ), ', which is outside grid'
          print *, 'Value reset to', sngl( edge )
          rtrunc( icmp ) = edge
        end if
        call cutoff( rtrunc( icmp ) )
      end if
c
c set potential offset
      Phimax = Phispl( rmax, zero )
c revise energy range
      call cutoff( rtrunc( icmp ) )
      dfcns( 3, icmp ) = Phimax
c reset velocity dispersion for lowered isothermal model (dfm = W_0)
      if( impar( icmp ) .eq. 2 )then
        sigmai2 =
     +    ( dfcns( 3, icmp ) - Phispl( zero, zero ) ) / dfcns( 1, icmp )
        print *, 'Dispersion now', sngl( sigmai2 )
      end if
c force re-initialization of any DF tables at this stage
      tothalo = distfn( Emine, 0.d0 )
c assign halo mass to grid
      call masscl( .false. )
      call switch( kgrid )
      call massdf( sphrho )
      call mascmb
      call switch( kgrid )
c sum halo mass
      if( s3d )then
        n = nr( jgrid )
        j = 0
        l = 4 * s3ntm * ( n - 1 ) + 1
        tothalo = s3dfld( 1, jgrid ) + s3dfld( l, jgrid ) * s3rad( n )
        l = j
      else
        tothalo = 0.
        l = 0
        do i = 1, mesh( jgrid )
          tothalo = tothalo + grdmss( i, 1 )
        end do
        jrad( 1 ) = 1
        krad( 1 ) = nr( jgrid )
      end if
c mass array was in program units
      tothalo = tothalo / ( lscale**3 * ts**2 )
      if( totdisc .gt. 0. )then
        print '( ''Halo mass and fraction now'', 2f16.7 )', tothalo,
     +                                   tothalo / ( totdisc + tothalo )
      else
        print *, 'Halo mass now', tothalo
      end if
c adjust normalization factor for DF
      m1 = m2
      m2 = tothalo
      icall = icall + 1
      norm = cmpmas( icmp ) * fmass( icmp ) / tothalo
      if( icall .eq. 1 )then
        norm = .99
      else if( icall .lt. 4 )then
        dfcns( 2, icmp ) = .9 * norm * dfcns( 2, icmp )
      else
c Newton's method
        m0 = cmpmas( icmp ) * fmass( icmp )
        dfcns( 2, icmp ) = ( m0 - m2 ) * ( b1 - b2 ) / ( m1 - m2 ) + b2
      end if
      b1 = b2
      b2 = dfcns( 2, icmp )
c rescale halo density distribution and combine with disc
      if( s3d )then
        if( hybrid )then
          do i = 1, mesh( jgrid )
            s3dmss( i, 1, 1 ) = norm * s3dmss( i, 1, jgrid )
            s3dmss( i, 1, 2 ) = mshd( i )
          end do
        else
          do i = 1, mesh( jgrid )
            s3dmss( i, 1, jgrid ) =
     +                          norm * s3dmss( i, 1, jgrid ) + mshd( i )
          end do
        end if
c reset force terms
        call s3dswp
        call switch( kgrid )
      else if( lgrd )then
        do i = 1, mesh( jgrid )
          grdmss( i, 1 ) = norm * grdmss( i, 1 ) + mshd( i )
        end do
      else
        call crash( 'DENSET', 'Unrecognized option' )
      end if
      return
      end
