      subroutine modfit( type )
c  Copyright (C) 2015, Jerry Sellwood
c
c Main driving routine for mode fitting
c
c It selects the type of fit, the subset of data to be fitted, the number
c   of modes and scaling factor.  If an acceptable solution is found, it
c   then offers to draw the fit and/or the eigenvectors and to save the
c   fit in a file
      use aarrays
      implicit none
c
c calling argument
      character*4 type
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      include 'inc/jscmmn.f'
c
c external
      logical gtlogl
c
c local variables
      integer i, ifail, imode, j, k, mact, nmodes
      logical found
      real end, start, sumsq
c
      if( freq( 1, 1 ) + freq( 2, 1 ) .eq. 0. )then
        do i = 1, 5
          freq( 1, i ) = .1 * real( 6 - i )
          freq( 2, i ) = .01 * real( 6 - i )
        end do
      end if
c read in Fourier harmonic
      call choosm( type, mact, i )
      if( mact .lt. 0 )return
c determine whether growth rate or pattern speed is to be held fixed
      ptonly = gtlogl( 'Force zero growth rates?' )
      gronly = .false.
      if( .not. ptonly )gronly = gtlogl( 'Force zero pattern speeds?' )
c start a new fitting procedure
      call readin( type )
    1 call select( mact, type )
    2 nmodes = 0
      do while ( ( nmodes .le. 0 ) .or. ( nmodes .gt. mmodes ) )
        call gtintg( 'Enter number of modes to be fitted', nmodes )
      end do
c choose exponent for rescaling data only if non-zero growth rates allowed
    3 if( ptonly )then
        apodc = 0
      else
        call gtreal( 'Enter exponential scale factor', apodc )
      end if
c check whether solution is on file
      call mbase( mact, nmodes, sumsq, .true., found )
      if( found )then
        if( gtlogl(
     +    'Existing solution found, do you just want to use it?' ) )then
          call alpbet( nmodes, sumsq )
          ifail = -1
          go to 4
        end if
      end if
c revise initial guess if requested
      j = 1
      k = 2
      if( ptonly )k = 1
      if( gronly )j = 2
      do imode = 1, nmodes
        print '( '' Current first guess is'', 2f8.3 )',
     +                                    ( freq( i, imode ), i = j, k )
        if( gtlogl( 'Do you want to change it?' ) )then
          if( ptonly )then
            call gtreal( 'Enter frequency', freq( 1, imode ) )
          else if( gronly )then
            call gtreal( 'Enter growth rate', freq( 2, imode ) )
          else
            call gtreals( 'Enter frequencies', freq( 1, imode ), 2 )
          end if
        end if
      end do
c call fitting routine
      call fitmod( nmodes, kp - jp + 1, sumsq, ifail )
      if( ifail .ge. 0 )call alpbet( nmodes, sumsq )
c draw fit and modes
    4 if( .not. null )then
        if( gtlogl( 'Do you want to see the fit' ) )call drwfit(
     +                                             mact, nmodes, sumsq )
        start = tme( jt )
        end = tme( kt )
        if( gtlogl( 'Do you want to see the mode shapes?' ) )then
          call drmode( mact, nmodes, start, end )
          if( gtlogl( 'Do you want to save the mode shape(s)?' )
     +                               )call absave( mact, nmodes, sumsq )
        end if
      end if
c fitting routine failed
      if( ifail .ge. 0 )then
c save new solution
        if( gtlogl( 'Do you want to save this in mbase?' )
     +                )call mbase( mact, nmodes, sumsq, .false., found )
      end if
c read instruction for what to do next
      do while ( .true. )
        call gtintg( 'Enter 0 to stop,' //
     + ' 1 to try a new time-base, 2 to change number of modes' //
     + ' or 3 to try a new apodisation const', i )
        if( i .le. 0 )return
        if( i .le. 3 )go to ( 1, 2, 3 ), i
      end do
      end
