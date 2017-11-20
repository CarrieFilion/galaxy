      real*8 function rhoisp( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c returns the density of a spherical isotropic model at the requested radius
c   by first looking up the relative potential at that radius and then using
c   the routine to evaluated the density as a function of psi
c
c calling argument
      real*8 r
c
c externals
      real*8 isopsi, isrhop
c
c local variable
      real*8 psi
c
      psi = isopsi( r )
      rhoisp = isrhop( psi )
      return
      end
