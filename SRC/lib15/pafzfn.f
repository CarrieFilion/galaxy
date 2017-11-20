      real*8 function pafzfn( rs )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c Vertical force at the field point specified in the common block arising
c   from a constant surface density annular mass element located at the
c   point rs in the z = 0 plane
c
c calling argument
      real*8 rs
c
c common block
c
      common / feildp / rf, zf
      real*8 rf, zf
c
c local variables
      real*8 fr, fz, pot
c
      call rngfrc( rs, rf, zf, fr, fz, pot )
      pafzfn = fz * rs
      return
      end
