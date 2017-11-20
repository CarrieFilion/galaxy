      subroutine sphbja( jst, coords, nact, jlen, besc )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c computes contribution from current group of particles to sph. Bessel coeffs
c
c calling arguments
      integer jlen, jst, nact
      real besc( nact, 2, jlen ), coords( 6, jst )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
      include 'inc/sphrbz.f'
c
c external
      real*8 sphbsj
c
c local allocatable arrays
      real, allocatable :: bes(:,:)
      real, allocatable :: cmp(:), smp(:)
c
c local variables
      integer is, ipp, j, k, l, m, n
      real arg, cp, rc, sp, splm
      real*8 rs
c
      if( uqmass )call crash( 'SPHBJA',
     +                        'Unequal particle masses not programmed' )
c allocate space
      allocate ( bes( jnmax, jlmax + 1 ) )
      allocate ( cmp( jlmax + 1 ) )
      allocate ( smp( jlmax + 1 ) )
c work through group
      do is = 1, jst
c skip non-halo particles and particles outside rbess
        ipp = iflag( is )
        if( ( .not. disc( ipp ) ) .and.
     +      ( rr( is ) .lt. rbess ) )then
          npbess( ipp ) = npbess( ipp ) + 1
c compute Bessel functions - l is index of Bessel fn, n is radial order
          do l = 0, jlmax
            do n = 1, jnmax
              rs = zeros( n, l + 1 ) * rr( is ) / rbess
              bes( n, l + 1 ) = sphbsj( rs, l )
            end do
          end do
c make table of cos(-m\phi) and sin(-m\phi)
          rc = sqrt( coords( 1, is )**2 + coords( 2, is )**2 ) + 1.e-8
          cp = coords( 1, is ) / rc
          sp = coords( 2, is ) / rc
          cmp( 1 ) = 1
          smp( 1 ) = 0
          do m = 1, jlmax
            cmp( m + 1 ) = cmp( m ) * cp + smp( m ) * sp
            smp( m + 1 ) = smp( m ) * cp - cmp( m ) * sp
          end do
c get Plms
          arg = coords( 3, is ) / rr( is )
          call s3tplm( dble( arg ), .false. )
          j = 0
          k = 0
          do l = 0, jlmax
            do m = 0, l
              k = k + 1
              splm = plm( k ) * pwt( is )
              do n = 1, jnmax
                j = j + 1
                besc( ipp, 1, j ) = besc( ipp, 1, j ) +
     +                             bes( n, l + 1 ) * splm * cmp( m + 1 )
                besc( ipp, 2, j ) = besc( ipp, 2, j ) +
     +                             bes( n, l + 1 ) * splm * smp( m + 1 )
              end do
            end do
          end do
        end if
      end do
c return local allocated space
      deallocate ( bes )
      deallocate ( cmp )
      deallocate ( smp )
      return
      end
