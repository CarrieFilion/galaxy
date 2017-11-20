      module interfac
!  Copyright (C) 2014, Richard James
      
         interface
            subroutine accumulate(x, par, perp, fac_par, fac_perp)        
               double precision, intent(inout)   :: x(:)	     	  
               double precision, intent(in)	 :: par(:), perp(:)  	  
               double precision, intent(in)	 :: fac_par, fac_perp
            end subroutine accumulate
         end interface
	 
	 interface
            subroutine boundary_pack
            end subroutine boundary_pack
	 end interface
	 
	 interface
            subroutine boundary_potential
            end subroutine boundary_potential
         end interface

	 interface
            subroutine boundary_save(preserve_hor, preserve_sn,           &
&	                             preserve_we)
               use mesh_specification
	       double precision, intent(out) :: preserve_hor(n23,2),	  &   
&           					preserve_sn(n13,2),	  &   
&           					preserve_we(n12,2)
	    end subroutine boundary_save
	 end interface
	 
	 interface
            subroutine boundary_smooth
	    end subroutine boundary_smooth
	 end interface
	 
         interface
            subroutine check_all
            end subroutine check_all
         end interface

         interface
            subroutine check_all_potentials
            end subroutine check_all_potentials
         end interface

         interface
            subroutine check_all_sources
            end subroutine check_all_sources
         end interface

         interface
            subroutine check_green(mesh, inter)
               double precision, intent(in)  :: mesh(:)
               logical, intent(in)           :: inter
            end subroutine check_green
         end interface

	 interface
	    subroutine check_internal(mesh, density, inter, deviation)
	       double precision, intent(inout) :: mesh(:)
	       double precision, intent(in)    :: density(:)
	       double precision, intent(out)   :: deviation
	       logical, intent(in)             :: inter
	    end subroutine check_internal
	 end interface
	 
         interface
	    subroutine check_poisson	    
	    end subroutine check_poisson   
         end interface
	 
	 interface
	    subroutine coerce(x, shift)
	       integer, intent(in)             :: shift
	       double precision, intent(inout) :: x(:)
	    end subroutine coerce
	 end interface
	 
	 interface
            subroutine edge_parallel(x, y, m, n)
               integer, intent(in)	     :: m, n
               double precision, intent(in)  :: x(m*n,2)
               double precision, intent(out) :: y(n, 2)
            end subroutine edge_parallel
	 end interface
	 
         interface
            subroutine edge_perpendicular(x, y, m, n)
               integer, intent(in)           :: m, n
               double precision, intent(in)  :: x(m*n,2)
               double precision, intent(out) :: y(m, 2)
            end subroutine edge_perpendicular
         end interface

	 interface
	    subroutine fftcsa(n, no, group, w, unity, halve)
	       integer, intent(in)             :: n, no, group
	       logical, intent(in)             :: unity, halve
	       double precision, intent(inout) :: w(n*no)
	    end subroutine fftcsa
         end interface
	 
	 interface
	    subroutine fftcss(n, no, group, w, unity, halved)
	       integer, intent(in)  :: n, no, group
	       double precision     :: w(n*no)
	       logical, intent(in)  :: unity, halved
	    end subroutine fftcss
         end interface

         interface
	    subroutine fftset(n, m)
	       integer, intent(in)            :: n, m(:)
	    end subroutine fftset
         end interface

	 interface
            subroutine fftsina_3d(w)
	       use mesh_specification
	       double precision, intent(inout), target :: w(n123)
            end subroutine fftsina_3d
	 end interface
	 
	 interface
	    subroutine fftsins_3d(w)
	       use mesh_specification
	       double precision, intent(inout), target :: w(n123)
	    end subroutine fftsins_3d
	 end interface
	 
	 interface
	    subroutine fftsna(n, no, group, w, unity)
	       integer, intent(in)             :: n, no, group
	       logical, intent(in)             :: unity
	       double precision, intent(inout) :: w(n*no)
	    end subroutine fftsna
         end interface

	 interface
	    subroutine greens_function
	    end subroutine greens_function
	 end interface
	 
	 interface
	    double precision function grenf3(i, j, k)
	       integer, intent(in)   :: i, j, k
	    end function grenf3
	 end interface
	 
	 interface
            subroutine initialise_solver
	    end subroutine initialise_solver
	 end interface
	 
	 interface
	    subroutine poiss
	    end subroutine poiss
	 end interface
	 
	 interface
            subroutine propagate_parallel(x, y, m, n)
               integer, intent(in)             :: m, n
               double precision, intent(inout) :: x(m*n)
               double precision, intent(in)    :: y(m,2)
            end subroutine propagate_parallel
         end interface

         interface
            subroutine propagate_perpendicular(x, y, m, n)
               integer, intent(in)             :: m, n
               double precision, intent(inout) :: x(m*n)
               double precision, intent(in)    :: y(n,2)
            end subroutine propagate_perpendicular
         end interface

	 interface
	    double precision function ran1(idum)
	       integer, intent(inout) :: idum
	    end function ran1
	 end interface
	 
	 interface
	    subroutine read_configuration
	    end subroutine read_configuration
	 end interface
	 
	 interface
	    subroutine setup_conversions
	    end subroutine setup_conversions
	 end interface
	 
	 interface
	    subroutine setup_green
	    end subroutine setup_green
	 end interface
	 
	 interface
	    subroutine setup_solver
	    end subroutine setup_solver
	 end interface
	 
	 interface
	    subroutine solve_laplace
	    end subroutine solve_laplace
	 end interface
	 
	 interface
	    subroutine solve_poisson
	    end subroutine solve_poisson
	 end interface
	 
	 interface
	    subroutine transp(x, nsets, ncont)
	       integer, intent(in) :: nsets, ncont
	       double precision    :: x(nsets*ncont)
	    end subroutine transp
	 end interface
	 
	 interface
	    subroutine vfthil(n, no, x, stc)
	       integer          :: n, no
	       logical          :: stc
	       double precision :: x(n*no)
	    end subroutine vfthil
	 end interface
	 
	 interface
	    subroutine vftsna(n, no, group, x)
	       integer, intent(in)             :: n, no, group
	       double precision, intent(inout) :: x(n*no)
	    end subroutine vftsna
	 end interface
	 
	 interface
	    subroutine vftsns(n, no, x)
	       integer, intent(in)             :: n, no
	       double precision, intent(inout) :: x(n*no)
	    end subroutine vftsns
	 end interface
	 
      end module interfac
