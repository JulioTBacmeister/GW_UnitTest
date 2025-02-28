module gw_interface
    use iso_c_binding
    implicit none
contains
    subroutine call_gw_rdg(arg1, arg2) bind(C, name="call_gw_rdg")
        use gw_common  ! Assuming gw_rdg depends on gw_common
        real(c_double), intent(in) :: arg1
        real(c_double), intent(out) :: arg2
        call gw_rdg(arg1, arg2)
    end subroutine call_gw_rdg
end module gw_interface
