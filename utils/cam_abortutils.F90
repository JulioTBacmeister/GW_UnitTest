module cam_abortutils

   !use shr_kind_mod, only: SHR_KIND_CL
   !use shr_sys_mod,  only: endrun => shr_sys_abort

   implicit none
   private
   save

   public :: endrun
   !public :: handle_allocate_error

CONTAINS

   !!subroutine handle_allocate_error(retval, subname, fieldname)
   subroutine endrun(textmess)
      ! if <retval> is not zero, generate an error message and abort
      ! Dummy arguments
      character(len=*), intent(in) :: textmess
      ! Local variable
      character(len=40)   :: errmsg

         write(errmsg, '(1a)') trim(textmess)
         write(*,*) errmsg
         STOP
         
    end subroutine endrun

end module cam_abortutils
