      subroutine choosm( type, mact, m )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to select sectoral harmonic for analysis from those available
c
c calling arguments
      character type*4
      integer m, mact
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/model.f'
c
c externals
      character*4 datnam
c
c local array
      integer maxk
      parameter ( maxk = 10 )
      integer sfpm( maxk )
c
c local variables
      character prompt*80
      integer i, ifail, j, k
      logical new, tryagain
c
c ensure logical is set
      do i = 1, ndattyp
        dattyp( i ) = dattyp( i ) .or. ( type .eq. datnam( i ) )
      end do
c get first header
      call rewnd
      call nextrec( type, ifail )
      if( ifail .ne. 0 )then
        print *,
     +     'Perhaps problem is population is not selected, icmp =', icmp
        call crash( 'CHOOSM', 'Data type not found' )
      end if
c determine available harmonics
      if( sfpc )then
        k = 1
        sfpm( 1 ) = msel( 1 )
        do i = 1, lastf
          new = .true.
          do j = 1, k
            if( msel( i ) .eq. sfpm( j ) )new = .false.
          end do
          if( new )then
            k = k + 1
            if( k .gt. maxk )call crash( 'CHOOSM',
     +                                          'sfpm array too small' )
            sfpm( k ) = msel( i )
          end if
        end do
      else
        k = 0
      end if
c choose sectoral harmonic
      tryagain = .true.
      do while ( tryagain )
        i = -10
        if( lgsp .or. sfpc .or. vlgs )i = 0
        if( danl .or. dnst .or. s3dc .or. scfc .or.
     +      sphb .or. vfld .or. zanl )i = -1
        if( i .lt. -1 )call crash( 'CHOOSM', 'Unrecognized data type' )
        j = ma * nsect
        if( danl )j = ( ma - 1 ) * nsect
        if( dnst .or. vfld )j = 0
        if( s3dc .or. scfc )j = s3lmax
        if( sfpc )then
          j = -1
          do i = 1, k
            j = max( j, sfpm( i ) )
          end do
          i = -1
        end if
        if( sphb )j = mr
        write( prompt, '( a, i3, a, i3, a )' )
     +        'Enter sectoral harmonic, up to', j, ' or', i, ' to stop'
        call gtintg( prompt, mact )
c force negative value to flag end
        if( mact .le. i )then
          mact = -1
          m = mact
          tryagain = .false.
        else
c ensure selected m is available
          tryagain = mact .gt. j
          if( sfpc )then
            m = mact
            tryagain = .true.
            do i = 1, k
              tryagain = tryagain .and. ( m .ne. sfpm( i ) )
            end do
          else
c all sectoral harmonics available for SPHB (and uncompressed ZANL data!)
            if( sphb )then
              m = mact
            else
c ensure selected m is a multiple of sectors
              m = mact / nsect
              i = m * nsect
              tryagain = tryagain .or. ( i .ne. mact )
            end if
          end if
        end if
      end do
      if( scfc )ncp = ( nbas + 1 ) * ( mact + 1 )
c (re-)select grid appropriate for population
      if( hybrid )then
        i = igrd( icmp )
        call switch( i )
      end if
      return
      end
