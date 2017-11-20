      character*4 function datnam( itype )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Returns a 4 character key word from one of the data types in the
c   runXXX.res file.  The input integer must be within the range
c   1 =< itype =< mtype
c
c calling argument
      integer itype
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
c local variables
      character*4 a( mtype )
      integer i, j, igt( 6 )
      logical firstc
c
      data firstc / .true. /
      data igt / 2, 11, 13, 19, 25, 28 /
      data a / 'PNTS', 'DNST', 'PNTL', 'FRQS', 'LGSP', 'ANGM', 'VFLD',
     +         'INTG', 'ORBS', 'GREY', 'DANL', 'MONI', 'SFPC', 'ZANL',
     +         'VFLH', 'SPHB', 'ZPRF', 'RNGS', 'SCFC', 'RHOR', 'SATE',
     +         'VLGS', 'MOIT', 'SIGR', 'S3DC', 'XCEN', 'LVAL', 'PROJ' /
c     +         'MODE', 'GASH', 'GPRP', 'PANL', 'CSPI'
c
c    1  'PNTS',    2  'DNST',    3  'PNTL',    4  'FRQS',    5  'LGSP',
c    6  'ANGM',    7  'VFLD',    8  'INTG',    9  'ORBS',   10  'GREY',
c   11  'DANL',   12  'MONI',   13  'SFPC',   14  'ZANL',   15  'VFLH',
c   16  'SPHB',   17  'ZPRF',   18  'RNGS',   19  'SCFC',   20  'RHOR',
c   21  'SATE',   22  'VLGS',   23  'MOIT',   24  'SIGR',   24  'S3DC',
c   26  'XCEN',   27  'LVAL',   28  'PROJ'
c
      if( ( itype .le. 0 ) .or. ( itype .gt. mtype )
     +                      )call crash( 'DATNAM', 'Unknown data type' )
c initialize common variables on first call
      if( firstc ) then
        ntype = mtype
c set flags for data types that depend on the grid and not the mass component
        do j = 1, mtype
          lgtyp( j ) = .false.
          do i = 1, 6
            if( igt( i ) .eq. j )lgtyp( j ) = .true.
          end do
        end do
        firstc = .false.
      end if
c return selected character string
      datnam = a( itype )
      return
      end
