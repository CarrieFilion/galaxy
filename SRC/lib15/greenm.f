      subroutine greenm( direct )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c driver routine to create the data tables needed for the field determination
c
c Grid methods generally require a Greens function.  The routine passes on
c   the calling argument, which indicates whether or not a direct version is
c   required in addition to the default transformed one.
c
c SFP methods need tables of function values and their derivatives.
c
c No tables are required for the PM+SH method.
c
c calling argument
      logical direct
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c local variable
      integer i
c
      if( parallel )call crash( 'GREENM', 'Parallel version needed' )
      do i = 1, ngrid
        call switch( i )
        if( ncode .gt. 11 )call crash( 'GREENM', 'Unrecognized method' )
c grid methods
        if( p2d )call p2grnm( direct )
        if( c3d )call c3grnm
        if( p3d )call p3grnm( direct )
        if( p3a )call pagrnm( direct )
        if( c2d )call c2grnm( direct )
      end do
      call switch( 0 )
      return
      end
