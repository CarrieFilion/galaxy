      subroutine filset
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to set up filter for Fourier harmonics
c    When m=0 terms are filtered out, the analytic radial acceleration
c    will be applied to every particle
c
c Called as part of grdset - further changes can be made to the filter in
c    START which will be preserved in the headers of the .res and .dmp files
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/bdpps.f'
c
      include 'inc/grids.f'
c
      include 'inc/lunits.f'
c
c external
      integer nchars
c
c local array
      integer l
      parameter ( l = 10 )
      integer sk( l )
c
c local variables
      character line*80
      integer i, j, k, m, n
      logical flag
      save flag, sk
c
      data sk / l * -1 /
c
      if( ncode .gt. 11 )call crash( 'FILSET', 'Unrecognized method' )
c read sectors and fixrad instruction only for first code
      if( ncode .eq. igrid( 1 ) )then
c read in sectors instruction - irrelevant for trees, p3a and Cartesian grids
        if( .not.
     +    ( noslfg .or. p3a .or. c2d .or. c3d .or. dr3d .or. bht ) )then
          call getline( ni, line )
          if( line( 1:4 ) .ne. 'SECT' )then
            print *, 'Next card is'
            print *, line
            call crash( 'FILSET', 'nsect card not found' )
          end if
          read( line( 11:80 ), *, err = 1 )nsect
        else
          nsect = 1
        end if
c fixed axisymmetric force - impossible for Cartesian grids and tree codes
        if( c2d .or. c3d .or. dr3d .or. bht )then
          fixrad = .false.
c essential if no self-G
        else if( noslfg )then
          fixrad = .true.
        else
c set by minm for SFP methods
          if( .not. ( sf2d .or. sf3d ) )then
c read instruction for other methods
            call getline( ni, line )
            if( line( 1:4 ) .ne. 'SKIP' )then
              print *, 'Next card is'
              print *, line
              call crash( 'FILSET', 'skip card not found' )
            end if
            n = nchars( line( 11:80 ) )
c determine if there are any non-blank entries
            flag = n .gt. 0
            if( flag )then
c parse line
              j = 0
              i = 10
              n = n + 10
              do while ( i .lt. n )
c skip leading blanks
                i = i + 1
                if( line( i:i ) .ne. ' ' )then
c count non-blank characters
                  k = 0
    2             if( i + k .lt. n )then
                    if( line( i+k:i+k ) .ne. ' ' )then
                      k = k + 1
                      go to 2
                    end if
                  end if
c read next entry
                  j = j + 1
                  read( line( i:i+k ), *, err = 1 )sk( j )
                  i = i + k
                end if
              end do
            end if
          end if
        end if
      end if
c set active functions
      if( sf2d .or. sf3d )then
        lastf = 0
        fixrad = minm .gt. 0
        ng = maxm + 1
        do m = minm, maxm, nsect
          do n = minn, maxn
            lastf = lastf + 1
            if( lastf .gt. mostf )call crash( 'FILSET',
     +                                               'mostf too small' )
            msel( lastf ) = m
            nsel( lastf ) = n
          end do
        end do
      else
c set default values - note ng = maxm + 1
        if( p3a .or. c2d .or. c3d .or. dr3d .or. bht )ng = 1
        if( ( s3d .or. scf ) .and. ( ng .eq. 0 ) )ng = s3lmax + 1
        lg( 1 ) = .false.
        if( p2d .or. p3d )then
          do i = 2, na / 2 + 1
            lg( i ) = .true.
          end do
        end if
c azimuthal filter: true means sectoral harmonic is skipped (Note: m = i - 1)
        if( ng .gt. 1 )then
          do i = 2, ng
            if( nsect .eq. 1 )then
              lg( i ) = .false.
            else
              lg( i ) = mod( i - 1, nsect ) .ne. 0
            end if
          end do
        end if
      end if
c add extra filter flags, if any
      if( flag )then
        do i = 1, l
          if( sk( i ) .ge. 0 )lg( sk( i ) + 1 ) = .true.
        end do
        fixrad = lg( 1 )
      end if
      return
    1 print *, 'Problem with data card:'
      print *, line( 1:60 )
      call crash( 'FILSET', 'Error reading data' )
      stop
      end
