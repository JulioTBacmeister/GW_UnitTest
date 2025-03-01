module geopotential

!---------------------------------------------------------------------------------
! Compute geopotential from temperature or
! compute geopotential and temperature from dry static energy.
!
! The hydrostatic matrix elements must be consistent with the dynamics algorithm.
! The diagonal element is the itegration weight from interface k+1 to midpoint k.
! The offdiagonal element is the weight between interfaces.
!
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
! Hacked by: J. Bacmeister, Feb 2025 for use in GW Unit test
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public geopotential_t

contains

!===============================================================================
  subroutine geopotential_t( ncol , pver ,                   &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     )

!-----------------------------------------------------------------------
!
! Purpose:
! Compute the geopotential height (above the surface) at the midpoints and
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol                  ! Number of longitudes
    integer, intent(in) :: pver                  ! Number of layers

    real(r8), intent(in) :: piln (ncol,pver+1)  !   - Log interface pressures
    real(r8), intent(in) :: pmln (ncol,pver)    !   - Log midpoint pressures
    real(r8), intent(in) :: pint (ncol,pver+1)  !   - Interface pressures
    real(r8), intent(in) :: pmid (ncol,pver)    !   - Midpoint pressures
    real(r8), intent(in) :: pdel (ncol,pver)    !   - layer thickness
    real(r8), intent(in) :: rpdel(ncol,pver)    !   - inverse of layer thickness
    real(r8), intent(in) :: t    (ncol,pver)    !   - temperature
    !real(r8), intent(in) :: q    (:,:,:)  ! (pcols,pver,:)- tracers (moist mixing ratios)
    real(r8), intent(in) :: q    (ncol,pver)    !   - water vapor
    real(r8), intent(in) :: rair (ncol,pver)    !   - Gas constant for dry air
    real(r8), intent(in) :: gravit        !               - Acceleration of gravity
    real(r8), intent(in) :: zvir (ncol,pver)    !   - rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: zi(ncol,pver+1)    !   - Height above surface at interfaces
    real(r8), intent(out) :: zm(ncol,pver)      !   - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    logical  :: lagrang                 ! Lagrangian vertical coordinate flag
    integer  :: ixq                     ! state constituent array index for water vapor
    integer  :: pverp
    integer  :: i,k,idx                 ! Lon, level indices, water species index
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
    real(r8) :: qfac(ncol,pver)         ! factor to convert from wet to dry mixing ratio
    real(r8) :: sum_dry_mixing_ratio(ncol,pver)! sum of dry water mixing ratios
    

    !CCPP-required variables (not used):
    integer            :: errflg
    character(len=512) :: errmsg

!
!-----------------------------------------------------------------------
!
    !Determine index for water vapor mass mixing ratio ... NOPE (jtb)
    !call cnst_get_ind('Q', ixq)

    !
    ! original code for backwards compatability with FV and EUL
    ! (Hacked out jtb 2/2025)
    pverp = pver + 1
    
    !dry air gas constant over gravity
    rog(:ncol,:) = rair(:ncol,:) / gravit

    ! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

    ! Compute zi, zm from bottom up.
    ! Note, zi(i,k) is the interface above zm(i,k)
    do k = pver, 1, -1

       ! First set hydrostatic elements consistent with dynamics

       do i = 1,ncol
          hkl(i) = pdel(i,k) / pmid(i,k)
          hkk(i) = 0.5_r8 * hkl(i)
       end do

       ! Now compute tv, zm, zi

       do i = 1,ncol
          tvfac   = 1._r8 + zvir(i,k) * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
       end do
    end do

  end subroutine geopotential_t
end module geopotential
