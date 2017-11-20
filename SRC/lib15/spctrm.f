      subroutine spctrm
c  Copyright (C) 2015, Jerry Sellwood
c
c Determines and contours the power spectra of the time series of various
c   data types accumulated during the run and saved in the runXXX.res file
c
c Called from ANALYS
c Graphics routines are from JSPLOT
      use aarrays
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/grids.f'
c
c external
      logical gtlogl
c
c local variables
      character*4 b, type
      integer lt, m, mact, mp
      logical ok, power
c
c select type of data
      call datype( type )
c decide type of contour plot
      ok = .false.
      do while ( .not. ok )
        call gtchar( 'Contours of time series or power spectrum?', b )
        call uppercase( b )
        power = b .eq. 'POWE'
        ok = power .or. ( b .eq. 'TIME' )
      end do
c select sectoral harmonic
      mact = 0
      call choosm( type, mact, m )
      do while ( mact .ge. 0 )
c select subset of data to use
        call readin( type )
        call select( mact, type )
c choose rescaling exponent
        ok = .false.
        do while ( .not. ok )
          call gtreal( 'Enter rescaling exponent', apodc )
          if( ( apodc * ( tme( kt ) - tme( jt ) ) .gt. 10. ) )then
            print *, 'this large value ', apodc, 'could cause problems'
            ok = .not. gtlogl( 'do you want to change it' )
          else
            ok = .true.
          end if
        end do
        call expscl
c find power spectrum
        if( power )call ftrans
c contour time series or power spectrum
        lt = kt - jt + 1
        mp = kp - jp + 1
        call contour( lt, mp, mact, power )
c read new sectoral harmonic
        call choosm( type, mact, m )
      end do
      return
      end
