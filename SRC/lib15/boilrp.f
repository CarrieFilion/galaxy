      subroutine boilrp( nagwarn )
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c routine to output boiler plate messages.  It is called only at the start
c   of execution of programs
c
c calling argument
      logical nagwarn
c
c  common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      real version
c
      if( master )then
        print *, 'Main program of the GALAXY package.' //
     +           '  Copyright (C) 2015, Jerry Sellwood'
        print *, 'This program comes with ABSOLUTELY NO WARRANTY.' //
     +           '  It is free software'
        print *, 'and you are welcome to' //
     +           ' redistribute it under certain conditions.'
        print *, 'See <http://www.gnu.org/licenses/> for details.'
c
        print '( '' GALAXY version'', f7.2 )', version( 0. )
        if( nagwarn )then
          if( lnag )then
            print *, 'This executable assumes the NAG library' //
     +               ' has been linked'
          else
            print *, 'This executable assumes that the NAG library' //
     +               ' is unavailable.'
            print *, 'It will execute until a call to a NAG routine' //
     +               ' is made'
          end if
        end if
      end if
      return
      end
