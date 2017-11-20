      subroutine smmass
c  Copyright (C) 2015, Jerry Sellwood
c
c Routine to create a smooth density distribution on the grid from the
c   analytic mass distribution.  Useful to avoid random fluctuations in
c   the particle distribution.
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      external dskrho, rhoxyz, sphrho
c
c local variables
      integer i, jcmp
      logical revers
c
      if( hybrid .and. ( ncmp .gt. 1 ) )then
        revers = .not. disc( 1 )
      else
        revers = .false.
      end if
c clear mass arrays
      mzone = nzones
      call masscl( .false. )
c smooth density is set in zone 1 only
      mzone = 1
c set up a smooth mass distribution on master node only
      if( master )then
        do jcmp = 1, ncmp
          if( revers )then
            icmp = ncmp + 1 - jcmp
          else
            icmp = jcmp
          end if
c set icomp for a compressed halo
          if( cmprssd )then
            icomp = 0
            do i = 1, ncomp
              if( iccmp( i ) .eq. icmp )icomp = i
            end do
          end if
c skip rigid components
          if( .not. rigidp( icmp ) )then
            i = igrd( icmp )
            call switch( i )
            if( disc( icmp ) )then
c set disk mass - no dfn possible
              if( c3d )then
                call massdf( dskrho )
              else if( s3d )then
                call mass3d
              else
                call massdi
              end if
c assign mass to the other grid also if a hybrid method
              if( hybrid )then
                if( ngrid .gt. 2 )call crash( 'SMMASS',
     +                                          'Option not available' )
                call switch( 2 )
                if( c3d )then
                  call massdf( dskrho )
                else if( s3d )then
                  call mass3d
                  do i = 1, mesh( jgrid )
                    s3dmss( i, 1, 2 ) = s3dmss( i, 1, 1 )
                  end do
                else
                  call massdi
                end if
              end if
            else
c set halo mass
              if( ctype( icmp ) .eq. 'UNIS' )then
                call massdf( rhoxyz )
              else
                call massdf( sphrho )
              end if
            end if
          end if
        end do
      end if
c combine masses
      mzone = nzones
      call mascmb
      return
      end
