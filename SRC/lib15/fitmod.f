      subroutine fitmod( nmodes, mp, sumsq, ifail )
c  Copyright (C) 2015, Jerry Sellwood
c
c part of mode fitting software
c
c routine to set up and call the NAG minimization routine and to report
c   and save the best fit frequencies
c   uses Burkardt's routine SUMSL
      use aarrays
      implicit none
c
c calling arguments
      integer nmodes, mp, ifail
      real sumsq
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
c externals
      external lsfun2, lsgrd2, vfparm
c
c local arrays
      integer mds2
      parameter ( mds2 = 2 * mmodes )
      real*8 d( mds2 ), omega( mds2 )
      integer, allocatable :: iw(:)
      real*8, allocatable :: w(:)
c
c local variables
      integer imode, io, k, liw, lw, nf
      integer nm2, uiparm
      logical firstc
      real*8 sumsq2, urparm
      save firstc, io, iw, liw, lw, w
      data firstc / .true. /
c
c check space
      if( nmodes .gt. mmodes )
     +                    call space( mmodes, nmodes, 'work', 'FITMOD' )
c set constants
      nm2 = 2 * nmodes
      if( ptonly .or. gronly )nm2 = nmodes
c put in current best guess
      k = 0
      do imode = 1, nmodes
        k = k + 1
        d( k ) = 1
        if( gronly )then
          omega( k ) = freq( 2, imode ) - apodc
        else
          omega( k ) = freq( 1, imode )
          if( .not. ptonly )then
            k = k + 1
            d( k ) = 1
            omega( k ) = freq( 2, imode ) - apodc
          end if
        end if
      end do
c reset flag
      iabort = 0
c values for NAG
      if( firstc )then
        liw = 60
        lw = 71 + mds2 * ( mds2 + 15 ) / 2
        allocate ( iw( liw ) )
        allocate ( w( lw ) )
        call new_unit( io )
        firstc = .false.
      end if
c initialize
      call deflt( 2, iw, liw, lw, w )
c flag this as a cold start
      iw( 1 ) = 12
c parameters to suppress all output
      do io = 19, 24
        iw( io ) = 0
      end do
      iw( 21 ) = io
      iw( 23 ) = -1
c least squares fitting routine
      call sumsl( nm2, d, omega, lsfun2, lsgrd2, iw, liw, lw, w,
     +                     uiparm, urparm, vfparm )
      call lsfun2( nm2, omega, nf, sumsq2, uiparm, urparm, vfparm )
      close( io, status = 'delete' )
      if( lprint )print *, 'return status for sumsl', iw( 1 )
      ifail = iw( 1 )
      sumsq = sumsq2 / denf
c save result
      k = 0
      do imode = 1, nmodes
        k = k + 1
        if( gronly )then
          freq( 1, imode ) = 0
          freq( 2, imode ) = omega( k ) + apodc
          print '( '' Mode'', i3, '' growth rate'', f10.4 )',
     +          imode, freq( 2, imode )
        else
          freq( 1, imode ) = omega( k )
          if( ptonly )then
            freq( 2, imode ) = 0
            print '( '' Mode'', i3, '' pattern speed'', f10.4 )',
     +            imode, freq( 1, imode )
          else
            k = k + 1
            freq( 2, imode ) = omega( k ) + apodc
            print '( '' Mode'', i3, '' pattern speed'', f10.4, ' //
     +            ' '' growth rate'', f10.4 )',
     +            imode, freq( 1, imode ), freq( 2, imode )
          end if
        end if
      end do
      return
      end

      subroutine vfparm
      print *, 'vfparm called!'
      return
      end
