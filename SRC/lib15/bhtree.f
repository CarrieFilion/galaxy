      subroutine bhtree
c  Copyright (C) 2014, Jerry Sellwood
      use aarrays
      implicit none
c builds the Barnes-Hut tree from the particle positions
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local arrays
      real cmax( 3 ), cmin( 3 )
c
c local variables
      integer i, ilt, ip, jt, kp, ksub, kt, lt, nt
      real cbound, xt, yt
      equivalence ( jt, xt ), ( kt, yt )
c
c initialise new tree
      call bhboun
      do i = 1, bhmt
        itup( i ) = 0
      end do
      jtmax = 1
      tsize = 0
      do i = 1, 3
        tsize = max( tsize, bound( 2, i ) - bound( 1, i ) )
      end do
      do i = 1, bhnlev
        if( bhtol .gt. 0. )then
          bhbox( i ) = tsize / 2**( i - 1 )
          bhbox( i ) = ( bhbox( i ) / bhtol )**2
        else
          bhbox( i ) = tsize**2
        end if
      end do
c tree pointer - n.b. ip in this loop is larger by 1 than the conventional inext
      ilt = nwpp - 3
c work over particles
      do ip = 1, lpf, nwpp
c cell boundaries
        do i = 1, 3
          cmin( i ) = bound( 1, i )
          cmax( i ) = bound( 1, i ) + tsize
        end do
c start at stem of tree
        jt = 1
    1   if( itup( jt ) .eq. 0 )then
c fill empty cell
          itup( jt ) = -ip
          ptcls( ip + ilt ) = xt
        else
          if( itup( jt ) .lt. 0 )then
c new subdivision required
            kp = -itup( jt )
            itup( jt ) = jtmax + 1
            jtmax = jtmax + 8
            if( jtmax .gt. bhmt )then
              call crash( 'BHTREE', 'insufficient space in tree array' )
            end if
c downward link
            i = ( jtmax - 1 ) / 8
            itdown( i ) = jt
c place pre-existing particle in appropriate sub-cell
            ksub = 1
            do i = 1, 3
              if( ptcls( kp + i - 1 ) - .5 * ( cmin( i ) + cmax( i ) )
     +                                                     .gt. 0. )then
                ksub = ksub + i
                if( i .eq. 3 )ksub = ksub + 1
              end if
            end do
            kt = jtmax - 8 + ksub
            itup( kt ) = -kp
            ptcls( kp + ilt ) = yt
          end if
c work down chain
          jt = itup( jt )
c find appropriate sub-cell
          do i = 1, 3
            cbound = .5 * ( cmin( i ) + cmax( i ) )
            if( ( ptcls( ip + i - 1 ) - cbound ) .gt. 0. )then
              jt = jt + i
              if( i .eq. 3 )jt = jt + 1
              cmin( i ) = cbound
            else
              cmax( i ) = cbound
            end if
          end do
          go to 1
        end if
      end do
c work through tree grouping masses
      do kt = 1, jtmax
        jt = jtmax + 1 - kt
        if( itup( jt ) .lt. 0 )then
c new particle - now ip is consistent with inext convention
          ip = -itup( jt ) - 1
          do i = 1, 3
            bhcom( i, jt ) = ptcls( ip + i )
          end do
          if( uqmass )then
            bhcom( 4, jt ) = ptcls( ip + ncoor + 1 ) * pmass
          else
            bhcom( 4, jt ) = pmass
          end if
        else if( itup( jt ) .gt. 0 )then
c centre of mass of grouping at next level up
          do i = 1, 4
            bhcom( i, jt ) = 0
          end do
          nt = itup( jt )
          do lt = nt, nt + 7
            if( itup( lt ) .ne. 0 )then
              bhcom( 4, jt ) = bhcom( 4, jt ) + bhcom( 4, lt )
              do i = 1, 3
                bhcom( i, jt ) = bhcom( i, jt ) +
     +                                   bhcom( i, lt ) * bhcom( 4, lt )
              end do
            end if
          end do
          do i = 1, 3
            bhcom( i, jt ) = bhcom( i, jt ) / bhcom( 4, jt )
          end do
        end if
      end do
      return
      end
