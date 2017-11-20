      subroutine header( lcomp )
c  Copyright (C) 2015, Jerry Sellwood
c
c routine called at start of all analysis software to read the header
c   record from the .res file
      implicit none
c
c calling argument
      logical lcomp
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/anlys.f'
c
      include 'inc/model.f'
c
c local variable
      integer i, j
c
c input run number
      call gtintg( 'Enter run number', irun )
c allow for comparison runs if requested
      if( lcomp )call cmprun( -1 )
c read header records
      in = -1
      call opnfil( in, 'res', 'unformatted', 'old', 'seq', i )
      if( i .ne. 0 )call crash( 'HEADER', '.res file not found' )
      call hedrec( in, .true. )
c set number of header records
      nhead = 1 + ncmp
      do j = 1, ngrid
        call switch( j )
        if( .not. ( noslfg .or. dr3d ) )nhead = nhead + 1
      end do
      call switch( 0 )
c separated .res files no longer an option
      separt = .false.
      scale_set = .false.
c rewind files and set counters
      call rewnd
      return
      end
