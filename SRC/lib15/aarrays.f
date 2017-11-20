      module aarrays
c Copyright (C) 2015, Jerry Sellwood
c
c allocatable arrays that are accessible globally
c
c particle data
      real, allocatable, save :: ptcls( : )
c
c masses in current zones
      real, allocatable, save :: grdmss( :, : )
      real, allocatable, save :: s3dmss( :, :, : )
      real*8, allocatable, save :: sfpmss( :, : )
c
c masses in previous zones
      real, allocatable, save :: grdmoz( :, : )
      real, allocatable, save :: s3dmoz( :, :, : )
      real*8, allocatable, save :: sfpmoz( :, : )
c
c acceleration components and potentials
      real, allocatable, save :: grdfld( :, : )
      real, allocatable, save :: s3dfld( :, : )
      real, allocatable, save :: c3dfld( :, : )
      real, allocatable, save :: sfpfld( :, : )
c
c tables for field methods
      real*8, allocatable, save :: sfprad( : )
      real*8, allocatable, save :: sfplst( :, : )
c
c tree architecture
      integer, allocatable, save :: itdown( : )
      integer, allocatable, save :: itup( : )
      real, allocatable, save :: bhbox( : )
      real, allocatable, save :: bhcom( :, : )
c
c direct particle
      real, allocatable, save :: drpt( :, : )
c
c positions, masses, and energies of most bound particles for centering
      real, allocatable, save :: bindp( :, :, : )
c
c arrays for analysis
      real, allocatable, save :: sdata( :, : )
      real, allocatable, save :: tme( : )
      real, allocatable, save :: wres( : )
c
c array for orbit integration in smooth
      real*8, allocatable, save :: orbtab( : )
c
c tables for adiabatic compression and splines
      real*8, allocatable, save :: arad( : ), frin( : ), Phin( : )
      real*8, allocatable, save :: gmn( :, : ), rhon( :, : )
      real*8, allocatable, save :: gmnc( :, : ), glamda( :, : )
      real*8, allocatable, save :: rhonc( :, : ), rlamda( :, : )
c same potential for both
      real*8, allocatable, save :: frinc( : ), flamda( : )
      real*8, allocatable, save :: Phinc( : ), plamda( : )
c
      end module aarrays

      include '../RAJ/mesh_data.f90'
