      subroutine switch( j )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c subroutine to switch method types for multiple grid codes.
c
c Does nothing if there is only one grid, else all existing code logicals
c   are reset at every call to this routine
c
c When the calling argument is a positive integer, the switches are set
c   for only the method selected by the calling argument
c When the calling argument is zero, switches are set to activate all the
c   the methods in use
c
c calling argument
      integer j
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/grids.f'
c
c local variable
      integer i
c
      if( ( j .lt. 0 ) .or. ( j .gt. ngrid ) )then
        print *, j
        call crash( 'SWITCH', 'nonsense calling argument' )
      end if
      if( ngrid .eq. 1 )then
c reset default jgrid
        jgrid = 1
      else
c switch between grids in hybrid or multi-grid run
        do i = 1, ncodes
          gtype( i ) = .false.
        end do
        if( j .gt. 0 )then
c activate selected grid
          ncode = igrid( j )
          if( ncode .eq. 0 )call crash( 'SWITCH',
     +                             'No self-G in a multi-grid method?' )
          gtype( ncode ) = .true.
          jgrid = j
        else
c activate all methods and assume grid gtype( ngrid ) is the largest
          do i = 1, ngrid
            gtype( igrid( i ) ) = .true.
          end do
          jgrid = ngrid
        end if
      end if
c set generic logicals
      lgrd = c2d .or. c3d .or. p2d .or. p3a .or. p3d .or. s3d
      lsfp = scf .or. sf2d .or. sf3d
      ldrct = bht .or. dr3d
      return
      end
