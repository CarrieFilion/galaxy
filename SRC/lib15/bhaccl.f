      subroutine bhaccl( jst )
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
c computes the acceleration components acting on the current group of
c   particles that are computed using the Barnes-Hut tree algorithm
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
c externals
      real softm, sftpot
c
c local arrays
      integer mlev
      parameter ( mlev = 100 )
      integer mon( mlev )
      real x( 3 ), d( 3 )
c
c local variables
      integer i, ilev, ip, is, jt, kt
      real afac, d2, ds
c
      do is = 1, jst
c locate the particle in the tree
        kt = ltree( is )
        ip = loc( is ) + 1
        nskip( is ) = .true.
        if( ip .gt. 0 )then
c regular massive particle
          do i = 1, 3
            x( i ) = bhcom( i, kt )
          end do
        else
c test particles are flagged with -ve loc
          do i = 1, 3
            x( i ) = oldc( i, is )
          end do
c test whether it shares a cell with a massive particle
          if( bhcom( 4, kt ) .ne. 0 )then
            d2 = 0
            do i = 1, 3
              d( i ) = bhcom( i, kt ) - x( i )
              d2 = d2 + d( i )**2
            end do
c need force from adjacent massive particle in the same cell
            afac = 0
            if( d2 .gt. 0. )then
              ds = sqrt( d2 ) / softl
              afac = bhcom( 4, kt )
              if( ds .lt. 2. )afac = afac * softm( ds )
              afac = afac / d2**1.5
            end if
            do i = 1, 3
              acc( i, is ) = acc( i, is ) + afac * d( i )
            end do
            if( phys )gpot( is ) =
     +                        gpot( is ) + bhcom( 4, kt ) * sftpot( d2 )
          end if
        end if
c start at stem of tree
        jt = 1
        ilev = 1
        mon( ilev ) = 7
c
    1   if( ( itup( jt ) .ne. 0 ) .and. ( itup( jt ) .ne. -ip ) )then
c compute distance to com of grouping
          d2 = 0
          do i = 1, 3
            d( i ) = bhcom( i, jt ) - x( i )
            d2 = d2 + d( i )**2
          end do
          if( ( itup( jt ) .ge. 0 ) .and.
     +        ( bhbox( ilev ) .gt. d2 ) )then
c resolve group of particles
            ilev = ilev + 1
            if( ilev .gt. mlev )then
              print *, 'too many levels of subdivision'
              call crash( 'BHACCL', 'insufficient mlev in accel' )
            end if
            mon( ilev ) = 0
            jt = itup( jt )
            go to 1
          end if
c acceptable grouping or single particle - apply softening rule
          afac = 0
          if( d2 .gt. 0. )then
            ds = sqrt( d2 ) / softl
            afac = bhcom( 4, jt )
            if( ds .lt. 2. )afac = afac * softm( ds )
            afac = afac / d2**1.5
          end if
          do i = 1, 3
            acc( i, is ) = acc( i, is ) + afac * d( i )
          end do
          if( phys )gpot( is ) =
     +                        gpot( is ) + bhcom( 4, jt ) * sftpot( d2 )
        end if
c move on to next pointer at this level
    2   mon( ilev ) = mon( ilev ) + 1
        if( mon( ilev ) .lt. 8 )then
          jt = jt + 1
          go to 1
        end if
c level finished
        ilev = ilev - 1
        if( ilev .ne. 0 )then
          jt = itdown( ( jt + 6 ) / 8 )
          go to 2
        end if
      end do
      return
      end
