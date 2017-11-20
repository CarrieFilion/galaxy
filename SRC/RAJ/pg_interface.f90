      module pg_interface
!  Copyright (C) 2014, Richard James

      interface
        subroutine pgask(flag)
          logical, intent(in)            :: flag
        end subroutine pgask
      end interface

      interface
        subroutine pgaxis(opt, x1, y1, x2, y2, v1, v2, step, nsub,      &
&                         dmajl, dmajr, fmin, disp, orient)
          character(len=*), intent(in)   :: opt
          real, intent(in)               :: x1, y1, x2, y2, v1, v2, step
          real, intent(in)               :: dmajl, dmajr, fmin, disp, orient
          integer, intent(in)            :: nsub
        end subroutine pgaxis
      end interface

      interface
        integer function pgband(mode, posn, xref, yref, x, y, ch)
          integer, intent(in)            :: mode, posn
          real, intent(in)               :: xref, yref
          real, intent(inout)            :: x, y
          character, intent(out)         :: ch
        end function pgband
      end interface

      interface
        subroutine pgbox(xopt, xtick, nxsub, yopt, ytick, nysub)
          character(len=*), intent(in)   :: xopt, yopt
          real, intent(in)               :: xtick, ytick
          integer, intent(in)            :: nxsub, nysub
        end subroutine pgbox
      end interface

      interface
        subroutine pgclos
        end subroutine pgclos
      end interface

      interface
        integer function pgcurs(x, y, ch)
          real, intent(inout)            :: x, y
          character, intent(out)         :: ch
        end function pgcurs
      end interface

      interface
        subroutine pgend
        end subroutine pgend
      end interface

      interface
        subroutine pgenv(xmin, xmax, ymin, ymax, just, axis)
          real, intent(in)               :: xmin, xmax, ymin, ymax
          integer, intent(in)            :: just, axis
        end subroutine pgenv
      end interface

      interface
        subroutine pgeras
        end subroutine pgeras
      end interface

      interface
        subroutine pglab(xlbl, ylbl, toplbl)
          character(len=*), intent(in)   :: xlbl, ylbl, toplbl
        end subroutine pglab
      end interface

      interface
        subroutine pgline(n, xpts, ypts)
          integer, intent(in)            :: n
          real, intent(in)               :: xpts(n), ypts(n)
        end subroutine pgline
      end interface

      interface
        subroutine pgmtxt(side, disp, coord, fjust, text)
          character(len=*), intent(in)   :: side, text
          real, intent(in)               :: disp, coord, fjust
        end subroutine pgmtxt
      end interface

      interface
        integer function pgopen(string)
          character(len=*), intent(in)   :: string
        end function pgopen
      end interface

      interface
        subroutine pgpage
        end subroutine pgpage
      end interface

      interface
        subroutine pgpanl(ix, iy)
          integer, intent(in)            :: ix, iy
        end subroutine pgpanl
      end interface

      interface
        subroutine pgpap(width, aspect)
          real, intent(in)               :: width, aspect
        end subroutine pgpap
      end interface

      interface
        subroutine pgpt(n, xpts, ypts, symbol)
          integer, intent(in)            :: n, symbol
          real, intent(in)               :: xpts(n), ypts(n)          
        end subroutine pgpt
      end interface
      
      interface
         subroutine pgptxt(x, y, angle, fjust, text)
	    real                         :: x, y, angle, fjust
	    character(len=*)             :: text
	 end subroutine pgptxt
      end interface

      interface
        subroutine pgqch(size)
          real, intent(out)              :: size
        end subroutine pgqch
      end interface
      
      interface
         subroutine pgqtxt(x, y, angle, fjust, text, xbox, ybox)
	    real                         :: x, y, angle, fjust
	    character(len=*)             :: text
	    real                         :: xbox(4), ybox(4)
	 end subroutine pgqtxt
      end interface

      interface
        subroutine pgqid(id)
          integer, intent(out)           :: id
        end subroutine pgqid
      end interface
      
      interface
         subroutine pgqcs(units, xch, ych)
	    integer                      :: units
	    real                         :: xch, ych
	 end subroutine pgqcs
      end interface
      
      interface
         subroutine pgqvsz(units, x1, x2, y1, y2)
	    integer, intent(in)          :: units
	    real, intent(out)            :: x1, x2, y1, y2
	 end subroutine pgqvsz
      end interface

      interface
        subroutine pgscf(font)
          integer, intent(in)            :: font
        end subroutine pgscf
      end interface

      interface
        subroutine pgsch(size)
          real, intent(in)               :: size
        end subroutine pgsch
      end interface

      interface
        subroutine pgsci(ci)
          integer, intent(in)            :: ci
        end subroutine pgsci
      end interface

      interface
        subroutine pgslct(i)
          integer, intent(in)            :: i
        end subroutine pgslct
      end interface

      interface
        subroutine pgsls(ls)
          integer, intent(in)            :: ls
        end subroutine pgsls
      end interface

      interface
        subroutine pgslw(lw)
          integer, intent(in)            :: lw
        end subroutine pgslw
      end interface

      interface
        subroutine pgsubp(nxsub, nysub)
          integer, intent(in)            :: nxsub, nysub
        end subroutine pgsubp        
      end interface

      interface
         subroutine pgtext(x, y, text)
           real, intent(in)              :: x, y
           character(len=*), intent(in)  :: text
         end subroutine pgtext
      end interface

      interface
        subroutine pgupdt
        end subroutine pgupdt
      end interface

      end module pg_interface
