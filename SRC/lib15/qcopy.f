      logical function qcopy( name )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c part of the management of res files
c interatively asks the user whether data of the given type in the
c   the input file should be copied to the output file
c
c calling argument
      real name
c
c local variable
      character*1 yn
c
      print '( '' Copy '', a4, ''? (y/n)'' )', name
      read '( a1 )', yn
      call lowercase( yn )
      qcopy = yn .eq. 'y'
      return
      end
