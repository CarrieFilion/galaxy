      subroutine selpop( ip )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to select a population of massive particles for analysis
c
c calling argument
      integer ip
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c external
      integer iscomp
c
c local variables
      integer i
      logical resel
c
c count active populations
      ip = 0
      do i = 1, ncmp
        if( nsp( i ) .gt. 0 )ip = ip + 1
      end do
      if( rngs )ip = ip - 1
c abort if none found
      if( ip .le. 0 )then
        call crash( 'SELPOP', 'No active population found' )
      else if( ip .eq. 1 )then
c only one active population
        if( ncmp .gt. 1 )then
          do i = 1, ncmp
            if( ( nsp( i ) .gt. 0 ) .and.
     +          ( cmpmas( i ) .gt. 0. ) )icmp = i
          end do
        else
          icmp = 1
        end if
      else
c select which if more than 1
        print *, ip, ' active populations'
        do i = 1, ncmp
          if( nsp( i ) .gt. 0 )print *, i, ': type = ', ctype( i ),
     +             ', disk =', disc( i ), ', # particles', nsp( i )
        end do
        i = 0
        resel = .true.
        do while ( resel )
          call gtintg( 'choose which to plot', i )
          if( ( i .gt. 0 ) .and. ( i .le. ncmp ) )then
            resel = nsp( i ) .eq. 0
          end if
        end do
        ip = i
c also set icmp and icomp (if needed)
        icmp = ip
        if( ( .not. disc( ip ) ) .and. cmprssd )icomp = iscomp( icmp )
      end if
      return
      end
