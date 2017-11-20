      subroutine datype( type )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c utility routine to input the data type for many different forms of analysis
c
c calling argument
      character*4 type
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/model.f'
c
c externals
      character*4 datnam
      logical gtlogl
c
c local variables
      integer i, j, n( mcmp )
      logical set
c
c select type of data
    1 set = .false.
      do while ( .not. set )
        call gtchar( 'Enter data type required', type )
        call uppercase( type )
        do i = 1, ndattyp
          dattyp( i ) = type .eq. datnam( i )
          set = set .or. dattyp( i )
        end do
        set = set .or. ( type( 1:3 ) .eq. 'END' )
      end do
c vlgs data comes in two forms
      if( vlgs )avlgs = gtlogl( 'Use azimuthal velocity data?' )
c select population
      icmp = 0
      if( lgsp .or. sphb .or. zanl )then
        j = 0
c warning: this will include a rings population, if present
        do i = 1, ncmp
          if( nsp( i ) .gt. 0 )then
            if( ( ( lgsp .or. zanl ) .and. disc( i ) )  .or.
     +          ( sphb .and. .not. disc( i ) ) )then
              j = j + 1
              n( j ) = i
            end if
          end if
        end do
        if( j .gt. 1 )then
          print '( 1x, a, 5i3 )',
     +                 'Available populations are', ( n( i ), i = 1, j )
          call gtintg( 'Select population', icmp )
        else if( j .lt. 1 )then
          print *, 'No suitable population for selected data type'
          go to 1
        else
          icmp = n( 1 )
        end if
      else if( rhor .or. sigr )then
        icmp = 1
      else
c data types that cannot be distinguished by component
        do i = 1, ncmp
          if( ( nsp( i ) .gt. 0 ) .and.
     +        ( .not. testp( i ) ) )then
            if( ( danl .or. sfpc ) .and. disc( i ) )icmp = i
            if( ( scfc .or. s3dc ) .and. ( .not. disc( i ) ) )icmp = i
          end if
        end do
      end if
      if( icmp .le. 0 )call crash( 'DATYPE', 'No population selected' )
c select grid appropriate for population
      if( hybrid )then
        i = igrd( icmp )
        call switch( i )
      end if
      return
      end
