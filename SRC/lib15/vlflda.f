      subroutine vlflda( jst, coords, nact, lprop, prop, nLbin, nLa )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c adds contribution of current group of particles to velocity field bins
c   and angular momentum distributions
c called from ANLGRP
c   the time-centred coordinates are passed as a calling argument
c
c calling arguments - number of particles in this group and time-centred coords
      integer jst, lprop, nact, nLbin
      integer nLa( nact, nLbin )
      real coords( 6, jst ), prop( nact, lprop )
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
c local variables
      integer ib, ipp, is, j
      logical mp
      real Lz, rc, t, vr, vth
      include 'inc/pi.f'
c
c work through group
      do is = 1, jst
        ipp = iflag( is )
        mp = ( ipp .le. ncmp ) .and. ( pwt( is ) .gt. 0. )
c compute cylindrical radius and z angular momentum
        if( threed )then
          rc = sqrt( coords( 1, is )**2 + coords( 2, is )**2 )
          Lz = coords( 5, is ) * coords( 1, is ) -
     +        coords( 4, is ) * coords( 2, is )
        else
          rc = rr( is )
          Lz = coords( 4, is ) * coords( 1, is ) -
     +        coords( 3, is ) * coords( 2, is )
        end if
c ensure rc is greater than zero
        rc = rc + 1.e-10
c add particle into appropriate angular momentum distribution
        if( angm .and. ( ipp .gt. 0 ) .and. mp )then
          if( disc( ipp ) )then
            j = pwt( is ) * Lz * Lzfac
          else
            j = .5 * pwt( is ) * Lz * Lzfac
            j = j + maxLz / 2
          end if
          j = max( j, -1 )
          j = min( j, maxLz )
          j = j + 2
          nLa( ipp, j ) = nLa( ipp, j ) + 1
        end if
c add contribution to disc velocity field
        if( vfld .and. disc( ipp ) .and. mp )then
          if( twod )then
            ib = ( rc - rinh( 1 ) ) * drfac
            if( ib .lt. ndrad )then
              t = atan2( coords( 2, is ), coords( 1, is ) )
              if( t .lt. 0. )t = t + 2. * pi
              j = thfac * t
              j = min( j, ndang - 1 )
              ib = ib * idskp + j * nprop + 1
              prop( ipp, ib ) = prop( ipp, ib ) + pwt( is )
              vth = Lz / rc
              vr = ( coords( 3, is ) * coords( 1, is ) +
     +               coords( 4, is ) * coords( 2, is ) ) / rc
              prop( ipp, ib + 1 ) = prop( ipp, ib + 1 ) +
     +                                                   pwt( is ) * vth
              prop( ipp, ib + 2 ) = prop( ipp, ib + 2 ) +
     +                                                pwt( is ) * vth**2
              prop( ipp, ib + 3 ) = prop( ipp, ib + 3 ) +
     +                                                    pwt( is ) * vr
              prop( ipp, ib + 4 ) = prop( ipp, ib + 4 ) +
     +                                                 pwt( is ) * vr**2
            end if
          else if( threed )then
            ib = rc * drfac
            if( ib .lt. ndrad )then
              t = atan2( coords( 2, is ), coords( 1, is ) )
              if( t .lt. 0. )t= t + 2. * pi
              j = thfac * t
              j = min( j, ndang - 1 )
              ib = ib * idskp + j * nprop + 1
              vth = Lz / rc
              vr = ( coords( 4, is ) * coords( 1, is ) +
     +               coords( 5, is ) * coords( 2, is ) ) / rc
              prop( ipp, ib ) = prop( ipp, ib ) + pwt( is )
              prop( ipp, ib + 1 ) = prop( ipp, ib + 1 ) +
     +                                                  pwt( is ) * vth
              prop( ipp, ib + 2 ) = prop( ipp, ib + 2 ) +
     +                                            pwt( is ) * vth * vth
              prop( ipp, ib + 3 ) = prop( ipp, ib + 3 ) +
     +                                                   pwt( is ) * vr
              prop( ipp, ib + 4 ) = prop( ipp, ib + 4 ) +
     +                                              pwt( is ) * vr * vr
              prop( ipp, ib + 5 ) = prop( ipp, ib + 5 ) +
     +                                       pwt( is ) * coords( 3, is )
              prop( ipp, ib + 6 ) = prop( ipp, ib + 6 ) +
     +                                    pwt( is ) * coords( 3, is )**2
              prop( ipp, ib + 7 ) = prop( ipp, ib + 7 ) +
     +                                       pwt( is ) * coords( 6, is )
              prop( ipp, ib + 8 ) = prop( ipp, ib + 8 ) +
     +                                    pwt( is ) * coords( 6, is )**2
            end if
          end if
        end if
c add contribution to halo velocity field
        if( vflh .and. ( .not. disc( ipp ) ) .and. mp )then
          ib = rc * drfac
          if( ib .lt. nhrbin )then
            j = nint( coords( 3, is ) * hzfac ) + nhzbin / 2
            if( ( j .ge. 0 ) .and. ( j .lt. nhzbin ) )then
              ib = ib * ihskp + j * nprop + 1
              vth = Lz / rc
              vr = ( coords( 4, is ) * coords( 1, is ) +
     +               coords( 5, is ) * coords( 2, is ) ) / rc
              prop( ipp, ib ) = prop( ipp, ib ) + pwt( is )
              prop( ipp, ib + 1 ) = prop( ipp, ib + 1 ) +
     +                                                   pwt( is ) * vth
              prop( ipp, ib + 2 ) = prop( ipp, ib + 2 ) +
     +                                             pwt( is ) * vth * vth
              prop( ipp, ib + 3 ) = prop( ipp, ib + 3 ) +
     +                                                    pwt( is ) * vr
              prop( ipp, ib + 4 ) = prop( ipp, ib + 4 ) +
     +                                               pwt( is ) * vr * vr
              prop( ipp, ib + 5 ) = prop( ipp, ib + 5 ) +
     +                                       pwt( is ) * coords( 3, is )
              prop( ipp, ib + 6 ) = prop( ipp, ib + 6 ) +
     +                                    pwt( is ) * coords( 3, is )**2
              prop( ipp, ib + 7 ) = prop( ipp, ib + 7 ) +
     +                                       pwt( is ) * coords( 6, is )
              prop( ipp, ib + 8 ) = prop( ipp, ib + 8 ) +
     +                                    pwt( is ) * coords( 6, is )**2
            end if
          end if
        end if
      end do
      return
      end
