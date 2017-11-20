      real*8 function dfNFW( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Eddington's formula for an isotropic DF in a spherical model
c   formula (4-140b) of BT (p 237)
c   uses QUADPAK routine DQAWSE
c
c uses cubic-spline interpolation in a table that is set up at the
c   first call.  The function for NFW asymptotes to a power law for
c   low energies in the cusp, so this is factored out of the spline
c
c calling argument
      real*8 E
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      external dfedfn !, rhopsi
      real*8 Phitot, quad_als, splint2 !, sderiv2
c
c local arrays
      integer ntab
      parameter ( ntab = 201 )
      real*8, allocatable :: c(:), dftb(:), Etab(:), lamda(:)
c
c local variables
      integer i, ier, iuse
      real*8 E1, Emine0, Eps, epsabs, epsrel, Psimax, rm
      real*8 nghalf, sl, tiny, zero
      parameter ( nghalf = -.5, zero = 0, tiny = 1.d-20 )
      include 'inc/pi.f'
      save c, Etab, dftb, iuse, lamda, Emine0
c
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
        allocate ( dftb( ntab ) )
        allocate ( Etab( ntab ) )
        if( master )print *, 'DFNFW: Building a table'
c place outer cutoff at a very large radius
        rm = rmax
        rmax = 100 * rscale( icmp )
        Phimax = Phitot( rmax )
        Psimax = Phimax - Phitot( zero )
        do i = 1, ntab - 1
c concentrate values in the table at low E
          E1 = dble( i ) / dble( 2 * ntab )
          E1 = ( E1 / ( 1. - E1 ) )**2
          E1 = E1 * ( Emaxe - Emine )
          Etab( i ) = E1 + Emine
          Eps = Phimax - Etab( i )
          Eps = min( Eps, Psimax )
c adaptive quadrature rule for sqrt-type singularity at the end of the range
          epsabs = 1.d-10
          epsrel = epsabs
          dfNFW = quad_als( dfedfn, zero, Eps, zero, nghalf, 1, epsabs,
     +                      epsrel, ier )
c ier = 2 means requested accuracy was not achieved
          if( ier .ne. 0 .and. ier .ne. 2 )then
            print *, 'ier =', ier, ' from QUAD_ALS'
c            call crash( 'DFEINA', 'QUADPACK error' )
          end if
c add boundary term
c          arg = sderiv2( zero, rhopsi, zero, Psimax ) / sqrt( Eps )
c          dfNFW = dfNFW + arg
c normalize
          dfNFW = dfNFW / ( sqrt( 8.d0 ) * pi**2 )
c factor out the power law dependence
          if( dfNFW .gt. 0.d0 )then
            dftb( i ) = dfNFW * E1**2.5
          else if( dfNFW .eq. 0.d0 )then
            dftb( i ) = -5.
            print *, 'dfNFW 0 for', i, E1
          else
            print *, 'DF -ve for E =', i, sngl( Eps ), sngl( dfNFW )
            call crash( 'DFNFW', 'DF negative' )
          end if
        end do
c fit an interpolating spline
        allocate ( c( ntab + 4 ) )
        allocate ( lamda( ntab + 4 ) )
        E1 = Etab( 2 )
        dfNFW =
     +       splint2( Etab, dftb, ntab - 1, E1, lamda, c, .true. )
        iuse = 1
c restore original rmax
        rmax = rm
        phimax = Phitot( rmax )
c save original Emine
        Emine0 = Emine
        if( master )print *, 'DFNFW: Table ready'
      end if
c look up value in table
      if( E .gt. Etab( ntab - 3 ) )then
c extrapolate as a power law
        sl = log10( dftb( ntab - 3 ) / dftb( ntab - 2 ) ) /
     +                    log10( Etab( ntab - 3 ) / Etab( ntab - 2 ) )
        dfNFW = dftb( ntab - 3 ) * ( E / Etab( ntab - 3 ) )**sl
      else if( E .lt. Etab( 1 ) )then
        dfNFW = dftb( 1 )
      else
        dfNFW =
     +       splint2( Etab, dftb, ntab - 1, E, lamda, c, .false. )
      end if
      E1 = max( E - Emine0, tiny )
      dfNFW = dfNFW * E1**(-2.5)
      return
      end
