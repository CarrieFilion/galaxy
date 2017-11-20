      subroutine disply( jst, coords, nbplot, dispy )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c adds selected particles to plot file when called from ANLGRP.  It is
c   also called from MEASURE to initialize and complete records
c
c calling arguments
      integer jst, nbplot
      real coords( 6, jst ), dispy( nbplot )
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
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
c local variables
      character*4 bstr
      integer i, ipp, is, j, jb( mcmp ), jump( mcmp ), k, l, nbuf
      integer ncount( mcmp ), npplot
      parameter ( npplot = 100000 )
      real ahol
      save jump, nbuf, ncount
c
      data bstr / 'PNTS' /
c
c setting up section
c
      if( jst .eq. 0 )then
        inplot = nbplot + 1
        is = 1
c apportion the buffer for the sizes of the separate populations
        do ipp = 1, ncmp
          if( nsp( ipp ) .gt. 0 )then
            jb( ipp ) =
     +        nint( real( npplot ) * real( nsp( ipp ) ) / real( nbod ) )
            jb( ipp ) = min( nbod, jb( ipp ) )
          end if
        end do
c ensure that the number of particles selected can fit into the buffer
        do while ( ncoor * inplot .gt. nbplot )
          inplot = 0
c set jump parameter to choose a suitable number of particles
          do ipp = 1, ncmp
            if( nsp( ipp ) .gt. 0 )then
              jump( ipp ) = is * nsp( ipp ) / jb( ipp )
              jump( ipp ) = max( jump( ipp ), 1 )
              ncount( ipp ) = jump( ipp ) - 1
c determine number of particles to be written
              inplot = inplot + ( nsp( ipp ) - 1 ) / jump( ipp ) + 1
            end if
          end do
          if( master )print *,
     +        'Number of particles to be written to plot file =', inplot
          if( is .gt. 100 )then
            if( master )then
              print *, 'Number of particles in each component',
     +                                     ( nsp( ipp ), ipp = 1, ncmp )
              print *, 'jump parameter for each component',
     +                                    ( jump( ipp ), ipp = 1, ncmp )
              print *, "ncoor, inplot, nbplot", ncoor, inplot, nbplot
            end if
            call crash( 'DISPLY',
     +                'Logic error or display buffer is far too small' )
          end if
          is = 2 * is
        end do
        nbuf = 0
c
c write out data
c
      else if( jst .lt. 0 )then
        if( parallel )call crash( 'DISPLY', 'mpi version needed' )
c write control record
        if( master )then
          lbplot = 2000
          read( bstr, '( a4 )' )ahol
          write( nphys )irun, ahol, istep, nbuf, lbplot
          l = 0
          do j = 1, nbuf, lbplot
            k = l + 1
            l = l + lbplot
            l = min( l, nbuf )
            write( nphys )( dispy( i ), i = k, l )
          end do
          if( lprint )write( no,
     + '( ''Number of particles written to plot file ='', i10 )' )inplot
        end if
c
c work through group
c
      else
        do is = 1, jst
c skip massless particles
          ipp = iflag( is )
          if( ( ipp .lt. ncmp ) .or. ( .not. rngs ) )then
c skip unwanted particles
            ncount( ipp ) = ncount( ipp ) + 1
            if( ncount( ipp ) .ge. jump( ipp ) )then
              ncount( ipp ) = 0
c save selected particles for plotting
              if( twod )then
                do i = 1, 4
                  dispy( nbuf + i ) = coords( i, is )
                end do
                nbuf = nbuf + 4
              else
                do i = 1, 6
                  dispy( nbuf + i ) = coords( i, is )
                end do
c flag particles by population
                dispy( nbuf + 3 ) = dispy( nbuf + 3 ) +
     +                                ( iflag( is ) - 1 ) * zshift
                nbuf = nbuf + 6
              end if
            end if
c check space
            if( nbuf .gt. nbplot )then
              print *, nbuf, nbplot
              call crash ( 'DISPLY', 'Plot buffer too small' )
            end if
          end if
        end do
      end if
      return
      end
