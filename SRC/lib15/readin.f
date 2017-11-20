      subroutine readin( type )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine to read down the runXXX.res file to determine the number of
c   functions and data times available for analysis.
c Allowed input types are 'LGSP', 'DANL', 'SFPC', 'SFP3', 'SPHB' & 'ZANL'
      use aarrays
      implicit none
c
c calling arguments
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
      include 'inc/anlsis.f'
c
      include 'inc/bdpps.f'
c
c local variable
      integer ifail
c
c restart file
      call rewnd
c read first record
      nt = 0
      call nextrec( type, ifail )
      do while ( ifail .eq. 0 )
c determine number of functions for first moment only
        if( nt .eq. 0 )then
          if( .not. scfc )ncp = 0
          if( danl .or. dnst .or. lgsp .or. s3dc .or.
     +        vfld .or. vlgs .or. zanl )ncp = mr
          if( sfpc )ncp = ma
          if( sphb )then
            jlmax = mr
            jnmax = ma
            ncp = jnmax * ( ( jlmax + 1 ) * ( jlmax + 2 ) ) / 2
            call sphbjz( 2 * jlmax - 1 )
          end if
          if( ncp .eq. 0 )call crash( 'READIN',
     +                                        'Unrecognized data type' )
        end if
c note any change to maxr
        if( sfpc )then
          if( wres( 1 ) .ne. minmaxr )print *, 'Change of maxr from',
     +  minmaxr, ' to', wres( 1 ), ' for it =', nt + 1, ' at time', time
        end if
c count this and get next record
        nt = nt + 1
        if( nt .gt. mtm )call space( mtm, nt, '/ lgsp /', 'READIN' )
        tme( nt ) = time
        call nextrec( type, ifail )
      end do
      if( nt .eq. 0 )call crash( 'READIN', 'No records found' )
      return
      end
