program test_centroid
  ! -----------------------------------------------------------------------
  ! Standalone unit test for weighted_centroid and vorticity_centroid_levels.
  ! No CAM dependencies -- r8 defined locally, subroutines included below.
  ! Compile and run:
  !   gfortran -o test_centroid test_centroid.f90 && ./test_centroid
  ! -----------------------------------------------------------------------
  implicit none

  integer, parameter :: r8 = selected_real_kind(12)

  ! --- grid ---------------------------------------------------------------
  integer, parameter :: ncol = 1
  integer, parameter :: pver = 10

  ! Heights in metres, surface to TOA (CAM ordering: k=1 is top, k=pver surface)
  ! We fill them TOA-first so index 1 = highest altitude.
  real(r8) :: zm(ncol, pver)
  real(r8) :: vorticity(ncol, pver)

  ! --- outputs ------------------------------------------------------------
  real(r8) :: z_centroid(ncol), z_centroid_top(ncol)
  integer  :: k_centroid(ncol), k_centroid_top(ncol)
  integer  :: steer_lev(ncol), launch_lev(ncol)
  real(r8) :: z_steer(ncol), z_launch(ncol)

  ! --- misc ---------------------------------------------------------------
  real(r8) :: z_min_pbl, z_max_pbl
  real(r8) :: w(ncol,pver)
  integer  :: k

  ! -----------------------------------------------------------------------
  ! Define a 10-level column: TOA at k=1 (20 km) down to surface at k=10 (0 km).
  ! Heights evenly spaced 20,17.8,15.6,...,0 km converted to metres.
  ! -----------------------------------------------------------------------
  do k = 1, pver
     zm(1, k) = (pver - k) * (20000._r8 / (pver - 1))
  end do

  write(*,'(a)') '--- zm profile (m) ---'
  do k = 1, pver
     write(*,'(a,i3,a,f10.1)') '  k=', k, '  zm=', zm(1,k)
  end do

  ! -----------------------------------------------------------------------
  ! Synthetic vorticity: Gaussian peak around z~10 km (k=5 in 20km column)
  ! so the true centroid is near 10 km.
  ! -----------------------------------------------------------------------
  do k = 1, pver
     vorticity(1,k) = exp( -0.5_r8 * ((zm(1,k) - 10000._r8)/3000._r8)**2 )
  end do

  write(*,'(/,a)') '--- abs(vorticity) profile ---'
  do k = 1, pver
     write(*,'(a,i3,a,f8.4)') '  k=', k, '  |zeta|=', vorticity(1,k)
  end do

  ! -----------------------------------------------------------------------
  ! TEST 1: weighted_centroid, no bounds
  ! Expected: centroid near 10000 m
  ! -----------------------------------------------------------------------
  w = abs(vorticity)
  call weighted_centroid(w, zm, ncol, pver, -huge(1._r8), huge(1._r8), &
                          z_centroid, k_centroid)
  write(*,'(/,a)') '=== TEST 1: full profile, no bounds ==='
  write(*,'(a,f10.1,a)') '  z_centroid = ', z_centroid(1), ' m  (expect ~10000)'
  write(*,'(a,i3)')       '  k_centroid = ', k_centroid(1)

  ! -----------------------------------------------------------------------
  ! TEST 2: weighted_centroid, z_min=2000 m (exclude surface/PBL)
  ! Expected: centroid still near 10000 m (Gaussian tails below 2 km are tiny)
  ! -----------------------------------------------------------------------
  call weighted_centroid(w, zm, ncol, pver, 2000._r8, huge(1._r8), &
                          z_centroid, k_centroid)
  write(*,'(/,a)') '=== TEST 2: z_min=2000 m (exclude PBL) ==='
  write(*,'(a,f10.1,a)') '  z_centroid = ', z_centroid(1), ' m  (expect ~10000)'
  write(*,'(a,i3)')       '  k_centroid = ', k_centroid(1)

  ! -----------------------------------------------------------------------
  ! TEST 3: weighted_centroid, z_max = 8000 m (exclude upper column)
  ! Expected: centroid pulled well below 10 km
  ! -----------------------------------------------------------------------
  call weighted_centroid(w, zm, ncol, pver, 2000._r8, 8000._r8, &
                          z_centroid, k_centroid)
  write(*,'(/,a)') '=== TEST 3: z_min=2000 m, z_max=8000 m ==='
  write(*,'(a,f10.1,a)') '  z_centroid = ', z_centroid(1), ' m  (expect < 8000)'
  write(*,'(a,i3)')       '  k_centroid = ', k_centroid(1)

  ! -----------------------------------------------------------------------
  ! TEST 4: all-zero vorticity -> sentinel return
  ! Expected: z_centroid = -huge, k_centroid = -1
  ! -----------------------------------------------------------------------
  w = 0._r8
  call weighted_centroid(w, zm, ncol, pver, -huge(1._r8), huge(1._r8), &
                          z_centroid, k_centroid)
  write(*,'(/,a)') '=== TEST 4: all-zero weight (sentinel check) ==='
  write(*,'(a,l)')  '  z_centroid == -huge: ', z_centroid(1) == -huge(1._r8)
  write(*,'(a,i3)') '  k_centroid          : ', k_centroid(1)

  ! -----------------------------------------------------------------------
  ! TEST 5: two-stage centroid via vorticity_centroid_levels
  ! Gaussian peak at 10 km with PBL cutoff at 2 km.
  ! Expected:
  !   steering_level near k=5 (z~10 km, first centroid of full profile)
  !   launch_level   at a higher k (lower z, i.e. lower in the top half)
  !   z_launch < z_steer (launch level is higher altitude)
  ! -----------------------------------------------------------------------
  vorticity(1,:) = [ (exp(-0.5_r8*((zm(1,k)-10000._r8)/3000._r8)**2), k=1,pver) ]
  z_min_pbl = 2000._r8
  z_max_pbl = huge(1._r8)

  call vorticity_centroid_levels(vorticity, zm, ncol, pver, &
                                  z_min_pbl, z_max_pbl, &
                                  steer_lev, launch_lev, &
                                  z_steer, z_launch)

  write(*,'(/,a)') '=== TEST 5: two-stage centroid ==='
  write(*,'(a,i3,a,f10.1,a)') '  steering_level = ', steer_lev(1),  &
       '   z_steer  = ', z_steer(1),  ' m'
  write(*,'(a,i3,a,f10.1,a)') '  launch_level   = ', launch_lev(1), &
       '   z_launch = ', z_launch(1), ' m'
  write(*,'(a,l)') '  z_launch > z_steer (launch above steering): ', &
       z_launch(1) > z_steer(1)

  write(*,'(/,a)') 'All tests complete.'

contains

  ! -----------------------------------------------------------------------
  ! Subroutines under test (copied verbatim from gw_movmtn_centroid_levels)
  ! -----------------------------------------------------------------------

  subroutine weighted_centroid(w, z, ncol, pver, z_min, z_max, z_centroid, k_centroid)
    integer,  intent(in)  :: ncol, pver
    real(r8), intent(in)  :: w(ncol,pver), z(ncol,pver)
    real(r8), intent(in)  :: z_min, z_max
    real(r8), intent(out) :: z_centroid(ncol)
    integer,  intent(out) :: k_centroid(ncol)
    real(r8) :: wi(pver), zi(pver), num, denom
    integer  :: i, k
    do i = 1, ncol
       zi = z(i,:)
       wi = w(i,:)
       where (zi < z_min .or. zi > z_max) wi = 0._r8
       num   = 0._r8
       denom = 0._r8
       do k = 1, pver-1
          num   = num   + 0.5_r8*(zi(k)*wi(k) + zi(k+1)*wi(k+1)) * (zi(k+1)-zi(k))
          denom = denom + 0.5_r8*(      wi(k) +         wi(k+1)) * (zi(k+1)-zi(k))
       end do
       if (denom /= 0._r8) then
          z_centroid(i) = num/denom
          k_centroid(i) = minloc(abs(zi - z_centroid(i)), dim=1)
       else
          z_centroid(i) = -huge(1._r8)
          k_centroid(i) = -1
       end if
    end do
  end subroutine weighted_centroid

  subroutine vorticity_centroid_levels(vorticity, zm, ncol, pver, &
       z_min_pbl, z_max_pbl, steering_level, launch_level, z_steer, z_launch)
    integer,  intent(in)  :: ncol, pver
    real(r8), intent(in)  :: vorticity(ncol,pver), zm(ncol,pver)
    real(r8), intent(in)  :: z_min_pbl, z_max_pbl
    integer,  intent(out) :: steering_level(ncol), launch_level(ncol)
    real(r8), intent(out) :: z_steer(ncol), z_launch(ncol)
    real(r8) :: absvort(ncol,pver), w_top(ncol,pver)
    real(r8) :: z_cent(ncol), z_top(ncol)
    integer  :: i
    absvort = abs(vorticity)
    call weighted_centroid(absvort, zm, ncol, pver, z_min_pbl, z_max_pbl, &
                            z_cent, steering_level)
    w_top = absvort
    do i = 1, ncol
       if (steering_level(i) > 0) then
          w_top(i, steering_level(i):pver) = 0._r8
       else
          w_top(i,:) = 0._r8
       end if
    end do
    call weighted_centroid(w_top, zm, ncol, pver, -huge(1._r8), huge(1._r8), &
                            z_top, launch_level)
    z_steer  = z_cent
    z_launch = z_top
  end subroutine vorticity_centroid_levels

end program test_centroid
