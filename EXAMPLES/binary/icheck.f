      subroutine icheck( jst, coords, lmoni, ehmon, llz, Lzval )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c Adds contributions from current group of particles to the global integrals
c Called from ANLGRP with the time-centred coordinates already computed and
c   passed as the second argument
c
c calling arguments
      integer jst, llz, lmoni
      real coords( 6, jst ), ehmon( 6, lmoni ), Lzval( llz )
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c externals
      real axipot, halpot, phcorr, satpot, tsfac
c
c local variables
      integer i, ipp, is, j
      real Lz!, E
      real ax, ay, ct, st, t, vx, vy, x, y
      real*8 sfac
c
c work through group
      do is = 1, jst
        ipp = iflag( is )
        if( ( ipp .gt. 0 ) .and. ( ipp .le. ncmp ) )then
          nspop( ipp ) = nspop( ipp ) + 1
          popm( ipp ) = popm( ipp ) + pwt( is )
c
          rr( is ) = 0
c          E = 0
          sfac = 1. / tsfac( iz( is ) )**2
c Cartesian coordinates
          do i = 1, ndimen
            rr( is ) = rr( is ) + coords( i, is )**2
            cm( i, ipp ) = cm( i, ipp ) + pwt( is ) * coords( i, is )
            pp( i, ipp ) =
     +               pp( i, ipp ) + pwt( is ) * coords( i + ndimen, is )
            for( i, ipp ) =
     +                   for( i, ipp ) + sfac * pwt( is ) * acc( i, is )
            claus = claus
     +               + coords( i, is ) * sfac * pwt( is ) * acc( i, is )
            ke( ipp ) =
     +               ke( ipp ) + pwt( is ) * coords( i + ndimen, is )**2
c            E = E + .5 * coords( i + ndimen, is )**2
          end do
          rr( is ) = sqrt( rr( is ) )
          pe( ipp ) = pe( ipp ) + .5 * pwt( is ) * gpot( is )
c          E = E + gpot( is )
          if( pertbn )pe( ipp ) = pe( ipp ) + pwt( is ) * satpot( is )
          if( fixrad )then
            pe( ipp ) = pe( ipp ) + pwt( is ) * axipot( rr( is ) )
c            E = E + axipot( rr( is ) )
          else
            if( suppl  )then
              pe( ipp ) = pe( ipp ) + pwt( is ) * phcorr( rr( is ) )
c              E = E + phcorr( rr( is ) )
            end if
            if( rigidh )then
              pe( ipp ) = pe( ipp ) + pwt( is ) * halpot( rr( is ) )
c              E = E + halpot( rr( is ) )
            end if
          end if
c angular momenta
          Lz = coords( ndimen + 2, is ) * coords( 1, is ) -
     +         coords( ndimen + 1, is ) * coords( 2, is )
          ang( 1, ipp ) = ang( 1, ipp ) + pwt( is ) * Lz
          if( threed )then
            ang( 2, ipp ) = ang( 2, ipp ) +
     +                 pwt( is ) * ( coords( 4, is ) * coords( 3, is ) -
     +                               coords( 6, is ) * coords( 1, is ) )
            ang( 3, ipp ) = ang( 3, ipp ) +
     +                 pwt( is ) * ( coords( 6, is ) * coords( 2, is ) -
     +                               coords( 5, is ) * coords( 3, is ) )
          end if

          if( twod )then
            x = oldc( 1, is )
            y = oldc( 2, is )
            ax = acc( 1, is )
            ay = acc( 2, is )
            vx = coords( 3, is )
            vy = coords( 4, is )
            write( no, 200 )x, y, vx, vy, ax, ay
            write( 33 )istep, is, x, y, vx, vy, ax, ay
c            write( no, 201 )( w( i, is ), i = 1, 4 )
          else
            write( no, 200 )( oldc( i, is ), i = 1, 6 )
  200       format( 6f12.6 )
c            write( no, 201 )( w( i, is ), i = 1, 8 )
c 201        format( 8f9.4 )
            write( no, 200 )( acc( i, is ), i = 1, 3 ), gpot( is )
          end if

c save specific angular momentum of every particle if requested
          if( lval )then
            i = numprocs * loc( is ) / nwpp + 1 + myid
            Lzval( i ) = Lz * gvfac / lscale
c            Lzval( 1, i ) = Lz * gvfac / lscale
c            Lzval( 2, i ) = E * gvfac**2
          end if
c store specific energy and angular momentum if needed
          if( moni )then
            j = loc( is ) / nmskip
            if( ( j * nmskip .eq. loc( is ) ) .and.
     +          ( j .lt. nmonit ) )then
              j = j + 1
              ehmon( 1, j ) = 0
              do i = ndimen + 1, ncoor
                ehmon( 1, j ) = ehmon( 1, j ) + coords( i, is )**2
              end do
              ehmon( 1, j ) = .5 * ehmon( 1, j ) + gpot( is )
              if( suppl )ehmon( 1, j ) =
     +                                ehmon( 1, j ) + phcorr( rr( is ) )
              if( rigidh )ehmon( 1, j ) =
     +                                ehmon( 1, j ) + halpot( rr( is ) )
              if( threed )then
                ehmon( 2, j ) = coords( 5, is ) * coords( 1, is )
     +                        - coords( 4, is ) * coords( 2, is )
                ehmon( 3, j ) = coords( 4, is ) * coords( 3, is )
     +                        - coords( 6, is ) * coords( 1, is )
                ehmon( 4, j ) = coords( 6, is ) * coords( 2, is )
     +                        - coords( 5, is ) * coords( 3, is )
              else
                ehmon( 2, j ) = coords( 4, is ) * coords( 1, is ) -
     +                          coords( 3, is ) * coords( 2, is )
                ehmon( 3, j ) = coords( 1, is )
                ehmon( 4, j ) = coords( 2, is )
                ehmon( 5, j ) = coords( 3, is )
                ehmon( 6, j ) = coords( 4, is )
              end if
            end if
          end if
        end if
      end do
      return
      end
