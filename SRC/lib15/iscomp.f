      integer function iscomp( i )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling argument
      integer i
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c local variable
      integer j
c
      if( ( i .lt. 1 ) .or. ( i .gt. ncmp ) )then
        print *, 'i =', i
        call crash( 'ISCOMP', 'Nonsense calling argument' )
      end if
c default is no compressed component
      iscomp = 0
      j = 1
c search for a matching compressed component
      do while ( j .le. ncmp )
        if( i .eq. iccmp( j ) )iscomp = j
        j = j + 1
      end do
      return
      end

