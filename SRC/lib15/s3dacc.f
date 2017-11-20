      subroutine s3dacc( jst )
c  Copyright (C) 2016, Jerry Sellwood
c
c routine to evaluate the potential and acceleration components acting on
c   the current group of particles from both the interior and exterior
c   masses using a spherical harmonic expansion on the grid shells
c In order to define an approximately continuous force field, the interior
c   (exterior) mass increases (decreases) linearly from the shell just interior
c   to the particle outwards to the next shell
      use aarrays
      implicit none
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c local allocatable arrays
      real, allocatable, save :: Q( :,: ), rG( :,: )
      real, allocatable, save :: cmp( : ), smp( : )
c
c local array
      real sf( 2, 4 )
c
c local variables
      integer i, is, j, k, l, lm, m, n, nrg
      logical firstc
      real arg, cp, f1, f2, f3, rad, r2, s, sp, t, tiny, u, x, y, z
      parameter ( tiny = 1.e-10 )
      save firstc
      data firstc / .true. /
c
c allocate space
      if( firstc )then
        allocate ( Q( 2, s3ntm ) )
        allocate ( rG( 2, s3ntm ) )
        allocate ( cmp( 0:s3maxl ) )
        allocate ( smp( 0:s3maxl ) )
        firstc = .false.
      end if
c compute shell numbers and weights for interpolation
      call weight( jst, .true. )
      nrg = nr( jgrid )
c work through particles
      do is = 1, jst
c initialise
        x = oldc( 1, is ) - xcen( 1, jgrid )
        y = oldc( 2, is ) - xcen( 2, jgrid )
        z = oldc( 3, is ) - xcen( 3, jgrid )
        rad = rr( is ) + tiny
        k = ncl( is )
        i = 4 * s3ntm * ( k - 1 ) - 4
        j = 4 * s3ntm * k - 4
        lm = 0
        do l = 0, s3lmax
          do m = 0, l
            i = i + 4
            j = j + 4
            lm = lm + 1
c create a sub-table of real and imaginary parts of A_lm and B_lm
            if( hybrid .and. ( label( is ) .eq. 2 ) )then
c combine coeffs from both grids if needed
              do n = 1, 2
                sf( n, 1 ) = s3dfld( i + n, 1 ) + s3dfld( i + n, 2 )
                sf( n, 2 ) = s3dfld( j + n, 1 ) + s3dfld( j + n, 2 )
                sf( n, 3 ) = s3dfld( i + n + 2, 1 ) +
     +                       s3dfld( i + n + 2, 2 )
                sf( n, 4 ) = s3dfld( j + n + 2, 1 ) +
     +                       s3dfld( j + n + 2, 2 )
              end do
            else
c single grid
              do n = 1, 2
                sf( n, 1 ) = s3dfld( i + n, jgrid )
                sf( n, 2 ) = s3dfld( j + n, jgrid )
                sf( n, 3 ) = s3dfld( i + n + 2, jgrid )
                sf( n, 4 ) = s3dfld( j + n + 2, jgrid )
              end do
            end if
            if( k .lt. nrg )then
c combine real and imaginary parts of interior and exterior masses
              do n = 1, 2
                Q( n, lm ) = wt( 1, is ) * sf( n, 1 ) +
     +                       wt( 2, is ) * sf( n, 2 ) +
     +                       wt( 1, is ) * sf( n, 3 ) +
     +                       wt( 2, is ) * sf( n, 4 )
                rG( n, lm ) = -( l + 1 ) * ( wt( 1, is ) * sf( n, 1 ) +
     +                                       wt( 2, is ) * sf( n, 2 ) )
     +                               + l * ( wt( 1, is ) * sf( n, 3 ) +
     +                                       wt( 2, is ) * sf( n, 4 ) )
              end do
            else
c particle beyond outer edge - all mass is interior and wt( 1, is ) = 0
              do n = 1, 2
                Q( n, lm ) = wt( 2, is ) * sf( n, 2 )
                rG( n, lm ) = -( l + 1 ) * wt( 2, is ) * sf( n, 2 )
              end do
            end if
          end do
        end do
c make table of cos(m\phi) and sin(m\phi)
        r2 = x**2 + y**2 + tiny
        t = sqrt( r2 )
        cp = x / t
        sp = y / t
        cmp( 0 ) = 1
        smp( 0 ) = 0
        do m = 1, s3lmax
          cmp( m ) = cmp( m - 1 ) * cp - smp( m - 1 ) * sp
          smp( m ) = smp( m - 1 ) * cp + cmp( m - 1 ) * sp
        end do
c get Plms and their derivatives - NB they include sqrt((l-m)!/(l+m)!)
        arg = z / rad
        call s3tplm( dble( arg ), .true. )
c sum over l
        lm = 0
        do l = 0, s3lmax
c apply Fourier filter
          if( lg( l + 1 ) )then
            lm = lm + l + 1
          else
c sum over m
            do m = 0, l
              lm = lm + 1
c product of e^{im\phi) with the summations over sources - real part only
              s = dplm( lm ) * ( Q( 1, lm ) * cmp( m ) -
     +                           Q( 2, lm ) * smp( m ) )
              t = 2 * m * plm( lm ) * ( Q( 1, lm ) * smp( m ) +
     +                                  Q( 2, lm ) * cmp( m ) )
              u = plm( lm ) * ( rG( 1, lm ) * cmp( m ) -
     +                          rG( 2, lm ) * smp( m ) )
              if( m .ne. 0 )then
                s = 2 * s
                u = 2 * u
              end if
c acceleration components
             f1 = -s * ( x * z / rad ) + t * ( y * rad**2 ) / r2 + x * u
             f2 = -s * ( y * z / rad ) - t * ( x * rad**2 ) / r2 + y * u
              f3 = s * ( r2 / rad ) + z * u
              acc( 1, is ) = acc( 1, is ) + f1 / rad**2
              acc( 2, is ) = acc( 2, is ) + f2 / rad**2
              acc( 3, is ) = acc( 3, is ) + f3 / rad**2
c potential
              if( phys )then
                s = Q( 1, lm ) * cmp( m ) - Q( 2, lm ) * smp( m )
                if( m .ne. 0 )s = 2 * s
                gpot( is ) = gpot( is ) - plm( lm ) * s
              end if
c end sum over m
            end do
c end sum over l
          end if
        end do
c add monopole term from the mass on the central grid point
c        if( .not. lg( 1 ) )then
c          gpot( is ) = gpot( is ) - sf( 1, 2 ) / rad
c          s = sf( j + 1, 2 ) / rad**3
c          acc( 1, is ) = acc( 1, is ) - x * s
c          acc( 2, is ) = acc( 2, is ) - y * s
c          acc( 3, is ) = acc( 3, is ) - z * s
c        end if
      end do
      return
      end
