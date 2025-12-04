module utils_mod

!
! This module 
!
use shr_kind_mod, only: r8 => shr_kind_r8
use interpolate_data, only: lininterp


implicit none
private


! Public interface
public :: make_pressures
public :: make_nc_filename
public :: leovy_alpha

!==========================================================================
contains
!==========================================================================

subroutine make_pressures( ncol, pver, ntim, hyai, hybi, hyam, hybm, PS, pint, pmid )

  integer, intent(in) :: ncol,pver,ntim
  real(r8) , allocatable, intent(in)  :: hyai(:), hybi(:), hyam(:), hybm(:)
  real(r8) , allocatable, intent(in)  :: PS(:,:)
  ! OUT ...
  real(r8) , allocatable, intent(out) :: pint(:,:,:), pmid(:,:,:)

  ! local
  integer icol,k,n

  allocate( pmid(ncol,pver,ntim) , pint(ncol, pver+1,ntim) )

  do n=1,ntim
     do k=1,pver
        do icol=1,ncol
           pmid(icol,k,n) = hyam(k)*100000. + hybm(k)*PS(icol,n)
        end do
     end do
  end do

  do n=1,ntim
     do k=1,pver+1
        do icol=1,ncol
           pint(icol,k,n) = hyai(k)*100000. + hybi(k)*PS(icol,n)
        end do
     end do
  end do

end subroutine make_pressures

!==========================================================================
subroutine leovy_alpha( pver, pref_edge, alpha )

  integer,  intent(in) :: pver
  real(r8), intent(in) :: pref_edge(:)

  
  ! Interpolated Newtonian cooling coefficients.
  real(r8) , allocatable, intent(out) :: alpha(:)


  integer :: k
  
  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer, parameter :: nalph = 71
  real(r8) :: alpha0(nalph) = [ &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.10133333_r8,  0.104_r8,       &
       0.108_r8,       0.112_r8,       0.116_r8,       0.12066667_r8,  &
       0.126_r8,       0.132_r8,       0.138_r8,       0.144_r8,       &
       0.15133333_r8,  0.16_r8,        0.17_r8,        0.18_r8,        &
       0.19_r8,        0.19933333_r8,  0.208_r8,       0.216_r8,       &
       0.224_r8,       0.232_r8,       0.23466667_r8,  0.232_r8,       &
       0.224_r8,       0.216_r8,       0.208_r8,       0.20133333_r8,  &
       0.196_r8,       0.192_r8,       0.188_r8,       0.184_r8,       &
       0.18266667_r8,  0.184_r8,       0.188_r8,       0.192_r8,       &
       0.196_r8,       0.19333333_r8,  0.184_r8,       0.168_r8,       &
       0.152_r8,       0.136_r8,       0.12133333_r8,  0.108_r8,       &
       0.096_r8,       0.084_r8,       0.072_r8,       0.061_r8,       &
       0.051_r8,       0.042_r8,       0.033_r8,       0.024_r8,       &
       0.017666667_r8, 0.014_r8,       0.013_r8,       0.012_r8,       &
       0.011_r8,       0.010333333_r8, 0.01_r8,        0.01_r8,        &
       0.01_r8,        0.01_r8,        0.01_r8                         &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(r8) :: palph(nalph) = [ &
       2.06115E-06_r8, 2.74280E-06_r8, 3.64988E-06_r8, 4.85694E-06_r8, &
       6.46319E-06_r8, 8.60065E-06_r8, 1.14450E-05_r8, 1.52300E-05_r8, &
       2.02667E-05_r8, 2.69692E-05_r8, 3.58882E-05_r8, 4.77568E-05_r8, &
       6.35507E-05_r8, 8.45676E-05_r8, 0.000112535_r8, 0.000149752_r8, &
       0.000199277_r8, 0.000265180_r8, 0.000352878_r8, 0.000469579_r8, &
       0.000624875_r8, 0.000831529_r8, 0.00110653_r8,  0.00147247_r8,  &
       0.00195943_r8,  0.00260744_r8,  0.00346975_r8,  0.00461724_r8,  &
       0.00614421_r8,  0.00817618_r8,  0.0108801_r8,   0.0144783_r8,   &
       0.0192665_r8,   0.0256382_r8,   0.0341170_r8,   0.0453999_r8,   &
       0.0604142_r8,   0.0803939_r8,   0.106981_r8,    0.142361_r8,    &
       0.189442_r8,    0.252093_r8,    0.335463_r8,    0.446404_r8,    &
       0.594036_r8,    0.790490_r8,    1.05192_r8,     1.39980_r8,     &
       1.86273_r8,     2.47875_r8,     3.29851_r8,     4.38936_r8,     &
       5.84098_r8,     7.77266_r8,     10.3432_r8,     13.7638_r8,     &
       18.3156_r8,     24.3728_r8,     32.4332_r8,     43.1593_r8,     &
       57.4326_r8,     76.4263_r8,     101.701_r8,     135.335_r8,     &
       180.092_r8,     239.651_r8,     318.907_r8,     424.373_r8,     &
       564.718_r8,     751.477_r8,     1000._r8                        &
       ]


  allocate( alpha(pver+1) )
  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa
  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._r8
     alpha0(k) = max(alpha0(k), 1.e-6_r8)
     palph(k) = palph(k)*1.e2_r8
  end do

  ! interpolate to current vertical grid and obtain alpha
  call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , pver+1)

end subroutine leovy_alpha

!-----------------------------------------------
function make_nc_filename(ncdata_root, year, month, day, secday) result(ncfile)
   implicit none
   !---- Inputs ----
   character(len=*), intent(in) :: ncdata_root   ! root path & prefix
   integer,          intent(in) :: year, month, day, secday

   !---- Output ----
   character(len=:), allocatable :: ncfile

   !---- Locals ----
   character(len=32) :: datestr
   character(len=:), allocatable :: tmp

   ! Build the date/time tag: 2004-06-15-00000
   write(datestr, '(I4.4,"-",I2.2,"-",I2.2,"-",I5.5)') year, month, day, secday

   ! Compose the full filename
   tmp = trim(ncdata_root)//'.'//trim(datestr)//'.nc'
   ncfile = tmp
end function make_nc_filename
!-----------------------------------------------

end module utils_mod
