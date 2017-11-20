      subroutine rewnd
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c part of the analysis software
c rewinds the input .res file and skips forward ove the header informatio
c   to reposition at the start of the data records
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/anlys.f'
c
      include 'inc/model.f'
c
c external
      character*4 datnam
c
c local variables
      character*4 string
      integer i, itype
c
c rewind main .res file
      rewind in
      if( ntype .eq. 0 )string = datnam( 1 )
c rewind and skip headers of separated files
      if( separt )then
        do itype = 1, ntype
          if( lunit( itype ) )then
            rewind itype + 50
            do i = 1, nhead
              read( itype + 50 )
            end do
          end if
        end do
      else
c skip header of main .res file
        do i = 1, nhead
          read( in )
        end do
      end if
c reset last step counters
      do itype = 1, ntype
        do i = 1, ncmp
          lstp( i, itype ) = -1
        end do
      end do
      return
      end
