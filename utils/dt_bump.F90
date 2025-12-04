module dt_bump
  implicit none
  private
  public :: bump_ymdh

contains

  pure logical function is_leap_gregorian(y) result(leap)
    integer, intent(in) :: y
    leap = (mod(y,4)==0 .and. (mod(y,100)/=0 .or. mod(y,400)==0))
  end function is_leap_gregorian

  pure integer function dim_of_month(y, m, use_leap) result(nd)
    integer, intent(in) :: y, m
    logical, intent(in) :: use_leap
    integer, parameter :: dim_tab(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    nd = dim_tab(m)
    if (m==2 .and. use_leap .and. is_leap_gregorian(y)) nd = 29
  end function dim_of_month

  ! floor-style div/mod for hours so remainder is always in [0,23]
  pure subroutine divmod24(total_hours, q_days, r_hours)
    integer, intent(in)  :: total_hours
    integer, intent(out) :: q_days, r_hours
    if (total_hours >= 0) then
       q_days =  total_hours / 24
       r_hours = mod(total_hours, 24)
    else
       ! floor division for negatives:
       q_days = - ((-total_hours + 23) / 24)
       r_hours = total_hours - 24*q_days   ! now 0 <= r_hours < 24
    end if
  end subroutine divmod24

  pure subroutine bump_days(y, m, d, delta_days, use_leap)
    integer, intent(inout) :: y, m, d
    integer, intent(in)    :: delta_days
    logical, intent(in)    :: use_leap
    integer :: days, nd

    days = delta_days
    if (days > 0) then
       do
          nd = dim_of_month(y, m, use_leap)
          if (d + days <= nd) then
             d = d + days
             exit
          else
             days = days - (nd - d + 1)
             d = 1
             m = m + 1
             if (m > 12) then
                m = 1
                y = y + 1
             end if
          end if
       end do
    else if (days < 0) then
       do
          if (d + days >= 1) then
             d = d + days
             exit
          else
             days = days + d
             m = m - 1
             if (m < 1) then
                m = 12
                y = y - 1
             end if
             d = dim_of_month(y, m, use_leap)
          end if
       end do
    end if
  end subroutine bump_days

  !------------------------------
  ! Bump (y0,m0,d0,h0) by dt_hours.
  ! If use_leap is absent => defaults to .true. (Gregorian leap years).
  !------------------------------
  pure subroutine bump_ymdh(y0, m0, d0, h0, dt_hours, y1, m1, d1, h1, use_leap)
    integer, intent(in)  :: y0, m0, d0, h0
    integer, intent(in)  :: dt_hours
    integer, intent(out) :: y1, m1, d1, h1
    logical, intent(in), optional :: use_leap
    logical :: leapcal
    integer :: q_days, r_hours, Htot

    leapcal = .true.; if (present(use_leap)) leapcal = use_leap

    ! start from inputs
    y1 = y0; m1 = m0; d1 = d0

    ! carry hours -> days, keeping remainder in [0,23] for any sign of dt_hours
    Htot = h0 + dt_hours
    call divmod24(Htot, q_days, r_hours)
    h1 = r_hours

    ! bump the date by q_days (can be positive or negative)
    call bump_days(y1, m1, d1, q_days, leapcal)
  end subroutine bump_ymdh

end module dt_bump
