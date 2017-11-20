      integer function imdcon(k)
c
      integer k
c
c  ***  return integer machine-dependent constants  ***
c
c     ***  k = 1 means return standard output unit number.   ***
c     ***  k = 2 means return alternate output unit number.  ***
c     ***  k = 3 means return  input unit number.            ***
c          (note -- k = 2, 3 are used only by test programs.)
c
c  +++  port version follows...
c
      external i1mach
      integer i1mach
      integer mdperm(3)
      data mdperm(1)/2/, mdperm(2)/4/, mdperm(3)/1/
      imdcon = i1mach(mdperm(k))
c  +++  end of port version  +++
c
c  +++  non-port version follows...
c     integer mdcon(3)
c     data mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/
c     imdcon = mdcon(k)
c  +++  end of non-port version  +++
c
 999  return
c  ***  last card of imdcon follows  ***
      end
