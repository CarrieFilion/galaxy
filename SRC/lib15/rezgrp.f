      subroutine rezgrp( jst )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c routine to change the time step for particles that crossed a zone boundary
c   during the previous step
c
c calling argument
      integer jst
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/buffer.f'
c
      include 'inc/grids.f'
c
c externals
      integer zoneno
      real tsfac
c
c local variables
      integer is, kz
      real tfac
c
c check whether rezoning is possible
      if( ( ( mzone .gt. 1 ) .or. ( nguard .gt. 0 ) ) .and.
     +    ( ilist .lt. nlists ) )then
c work through all particles
        do is = 1, jst
c skip particles in outer zones
          if( iz( is ) .le. mzone )then
            jgrid = label( is )
            kz = zoneno( is, .true. )
            kz = min( kz, mzone )
c additional test for outward movers
            if( kz .gt. iz( is ) )kz = zoneno( is, .false. )
c test whether this is different from old zone
            if( kz .ne. iz( is ) )then
c get time step factor
              tfac = tsfac( kz ) / tsfac( iz( is ) )
c change time step - accelerations already available
              call chstep( is, tfac )
c keep track of number in central hole and intensive care
              if( ( kz .le. 0 ) .and. ( iz( is ) .gt. 0 ) )then
                nshortm = nshortm + 1
                if( hole )ncenm( iflag( is ) ) =
     +                                          ncenm( iflag( is ) ) + 1
              end if
              if( ( kz .gt. 0 ) .and. ( iz( is ) .le. 0 ) )then
                nshortm = nshortm - 1
                if( hole )ncenm( iflag( is ) ) =
     +                                          ncenm( iflag( is ) ) - 1
              end if
              iz( is ) = kz
            end if
          end if
        end do
      end if
      return
      end
