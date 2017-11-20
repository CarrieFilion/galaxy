      real*8 function dfEina( E )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c isotropic DF for the Einasto halo computed by Eddington's inversion formula
c   (4-140b) of BT (p 237)
c   uses QUADPAK routine DQAWSE
c
c uses cubic-spline interpolation in a table that is set up at the
c   first call.  The log of the function is tabulated to cope with
c   a wide range of values
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
      real*8 Psimin, Psimax
      common / dfeddl / Psimin, Psimax
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
      real*8 E1, Eps, epsabs, epsrel, nghalf, zero!, arg
      parameter ( nghalf = -.5, zero = 0 )
      include 'inc/pi.f'
      save c, Etab, dftb, iuse, lamda
c
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
        allocate ( dftb( ntab ) )
        allocate ( Etab( ntab ) )
        if( master )print *, 'DFEINA: Building a table'
        Psimin = zero
        Psimax = Phimax - Phitot( zero )
        do i = 1, ntab - 1
c concentrate values in the table at low E
          E1 = dble( i - 1 ) / dble( 2 * ( ntab - 1 ) )
          E1 = ( E1 / ( 1. - E1 ) )**2
          E1 = Emine + E1 * ( Emaxe - Emine )
          Etab( i ) = E1
          Eps = Phimax - E1
          Eps = min( Eps, Psimax )
c          Eps = min( Eps, Phimax - Phitot( 1.d-12 ) )
c adaptive quadrature rule for sqrt-type singularity at the end of the range
          epsabs = 1.d-10
          epsrel = epsabs
          dfEina = quad_als( dfedfn, zero, Eps, zero, nghalf, 1, epsabs,
     +                       epsrel, ier )
c ier = 2 means requested accuracy was not achieved
          if( ier .ne. 0 .and. ier .ne. 2 )then
            print *, 'ier =', ier, ' from QUAD_ALS'
c            call crash( 'DFEINA', 'QUADPACK error' )
          end if
c add boundary term
c          arg = sderiv2( zero, rhopsi, zero, Psimax ) / sqrt( Eps )
c          dfEina = dfEina + arg
c normalize
          dfEina = dfEina / ( sqrt( 8.d0 ) * pi**2 )
c save as the log, for better interpolation
          if( dfEina .gt. 0.d0 )then
            dftb( i ) = log10( dfEina )
          else if( dfEina .eq. 0.d0 )then
            dftb( i ) = -5.
            print *, 'dfEina 0 for', i, E1
          else
            print *, 'DF -ve for E =', sngl( E1 ), sngl( dfEina )
            call crash( 'DFEINA', 'DF negative' )
          end if
        end do
        dftb( ntab ) = -5.
        iuse = 1
        if( master )print *, 'DFEINA: Table ready'
c fit an interpolating spline
        allocate ( c( ntab + 4 ) )
        allocate ( lamda( ntab + 4 ) )
        E1 = Etab( 1 )
        dfEina = splint2( Etab, dftb, ntab, E1, lamda, c, .true. )
      end if
c look up value in table
      if( ( E .lt. Etab( 1 ) ) .or. ( E .gt. Etab( ntab ) ) )then
        dfEina = -20
        if( E - Etab( ntab ) .lt. 1.d-5 )dfEina = dftb( ntab )
      else
        dfEina = splint2( Etab, dftb, ntab, E, lamda, c, .false. )
      end if
      dfEina = 10.d0**dfEina
      return
      end
