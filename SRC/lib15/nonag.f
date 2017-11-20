      subroutine nonag( from, named )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to provide a graceful exit from execution when an unavalable
c   NAG routine is required.  No longer needed from v15
c
c calling arguments
      character*(*) from, named
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      if( master )then
        write( no, '( 1x, a )' )'Substitute for NAG routine ' //
     +         named // ' needed in subprogram ' // from
        print *, 'Substitute for NAG routine ' // named //
     +           ' needed in subprogram ' // from
        print *
        print *, 'If the user has access to the NAG library, then' //
     +           ' the GALAXY package should'
        print *, 'be remade as described in sec 1.4 of the manual'
        print *
        print *, 'If the user does not have NAG, then an equvalent' //
     +           ' routine will be needed.'
        print *, 'The calling arguments and functionality of the' //
     +           ' routine to be substituted'
        print *, 'are described in the NAG on-line documentation:'
        print *,
     + 'http://www.nag.co.uk/numeric/fl/manual/html/FLlibrarymanual.asp'
      end if
      call crash( 'NONAG', 'Unavailable NAG routine called' )
      end
