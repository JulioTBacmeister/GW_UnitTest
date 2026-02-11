module nc_flexout_mod
   use netcdf
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private

   integer, save :: ncid        = -1
   integer, save :: dimid_col   = -1
   integer, save :: dimid_x     = -1
   integer, save :: dimid_y     = -1
   integer, save :: dimid_z     = -1
   integer, save :: dimid_ze    = -1
   integer, save :: dimid_time  = -1
   integer, save :: ncol_save   = 0
   integer, save :: nx_save     = 0
   integer, save :: ny_save     = 0
   integer, save :: nz_save     = 0
   integer, save :: nze_save    = 0


   public :: ncfile_init
   public :: ncfile_init_col
   public :: ncfile_set_globals
   public :: ncfile_put_3d
   public :: ncfile_put_col1d_notime
   public :: ncfile_put_col2d_notime
   public :: ncfile_put_col2d
   public :: ncfile_put_col3d
   public :: ncfile_close
   
contains

   !----------------------------
   subroutine check_nc(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where
      if (ierr /= nf90_noerr) then
         write(*,*) 'NetCDF error in ', trim(where), ': ', trim(nf90_strerror(ierr))
         stop 1
      end if
   end subroutine check_nc

   !----------------------------
   subroutine ncfile_init(filename, nx, ny, nz)
      character(len=*), intent(in) :: filename
      integer,          intent(in) :: nx, ny, nz
      integer :: ierr

      ! Create file (NETCDF4 so redef is cheap; use NF90_CLOBBER|NF90_NETCDF4 if desired)
      ierr = nf90_create(filename, NF90_NETCDF4, ncid)
      call check_nc(ierr, 'nf90_create')

      ! Define dims; time is unlimited
      ierr = nf90_def_dim(ncid, 'lon',   nx,          dimid_x);     call check_nc(ierr, 'def_dim lon')
      ierr = nf90_def_dim(ncid, 'lat',   ny,          dimid_y);     call check_nc(ierr, 'def_dim lat')
      ierr = nf90_def_dim(ncid, 'level', nz,          dimid_z);     call check_nc(ierr, 'def_dim level')
      ierr = nf90_def_dim(ncid, 'time',  NF90_UNLIMITED, dimid_time); call check_nc(ierr, 'def_dim time')

      ! (Optional) define coordinate variables here, attributes, etc.

      ierr = nf90_enddef(ncid)
      call check_nc(ierr, 'nf90_enddef(init)')

      nx_save = nx
      ny_save = ny
      nz_save = nz
   end subroutine ncfile_init

   !----------------------------
   subroutine ncfile_init_col(filename, ncol, nz, time_val, ny, nx, lat_R, lon_R )
      character(len=*),  intent(in) :: filename
      integer,           intent(in) :: ncol, nz
      real(r8),          intent(in) :: time_val
      integer,  optional, intent(in) :: ny, nx
      real(r8), optional, intent(in) :: lat_R(:), lon_R(:)
      integer :: ierr
      integer :: varid_lon_R, varid_lat_R, varid_time

      ! Create file (NETCDF4 so redef is cheap; use NF90_CLOBBER|NF90_NETCDF4 if desired)
      ierr = nf90_create(filename, NF90_NETCDF4, ncid)
      call check_nc(ierr, 'nf90_create')

      ! Define dims; time is unlimited
      ierr = nf90_def_dim(ncid, 'ncol',   ncol,        dimid_col);   call check_nc(ierr, 'def_dim col')
      ierr = nf90_def_dim(ncid, 'level', nz,          dimid_z);     call check_nc(ierr, 'def_dim level')
      ierr = nf90_def_dim(ncid, 'ilevel', nz+1,       dimid_ze);     call check_nc(ierr, 'def_dim ilevel')
      ierr = nf90_def_dim(ncid, 'time',  NF90_UNLIMITED, dimid_time); call check_nc(ierr, 'def_dim time')

      ! add time var - should appear as coordinate since name is same as dim
      ierr = nf90_def_var(ncid, 'time', NF90_DOUBLE, (/ dimid_time /), varid_time)

      if (present(ny)) then
         write(*,*) "NY is here",ny
         ierr = nf90_def_dim(ncid, 'ny',   ny,   dimid_y);   call check_nc(ierr, 'def_dim y')
      end if
      if (present(ny) .and.present(lat_R) ) then
         ierr = nf90_def_var(ncid, 'lat_R', NF90_DOUBLE, (/ dimid_y /), varid_lat_R)
         call check_nc(ierr, 'def_var lat_R')
      end if


      if (present(nx)) then
         write(*,*) "NX is here",nx
         ierr = nf90_def_dim(ncid, 'nx',   nx,   dimid_x);   call check_nc(ierr, 'def_dim x')
      end if
      if (present(nx) .and.present(lon_R) ) then
         ierr = nf90_def_var(ncid, 'lon_R', NF90_DOUBLE, (/ dimid_x /), varid_lon_R)
         call check_nc(ierr, 'def_var lon_R')
      end if

      ! (Optional) define coordinate variables here, attributes, etc.

      ierr = nf90_enddef(ncid)
      call check_nc(ierr, 'nf90_enddef(init)')

      
      ierr = nf90_put_var(ncid, varid_time, time_val)  !
      call check_nc(ierr, 'put_var time')
      
      if (present(ny) .and. present(lat_R)) then
         ierr = nf90_put_var(ncid, varid_lat_R, lat_R)  !
         call check_nc(ierr, 'put_var lat_R')
      end if
      if (present(nx) .and. present(lon_R)) then
         ierr = nf90_put_var(ncid, varid_lon_R, lon_R)  !
         call check_nc(ierr, 'put_var lon_R')
      end if

      ncol_save = ncol
      nz_save = nz
      nze_save = nz+1
    end subroutine ncfile_init_col
    
    !============================================================
    subroutine ncfile_set_globals(title, institution, source, &
                                  references, comment, ncdata, &
                                  calculation_type )
      !===========================
      !integer,          intent(in) :: ncid ! In preamble ... 
      character(*),     intent(in), optional :: title, institution, source, references, comment, ncdata, calculation_type
      integer :: ierr

      ! You must be in define mode to add/modify attributes
      ierr = nf90_redef(ncid)

      if (present(title))      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'title',       trim(title))
      if (present(institution))ierr = nf90_put_att(ncid, NF90_GLOBAL, 'institution', trim(institution))
      if (present(source))     ierr = nf90_put_att(ncid, NF90_GLOBAL, 'source',      trim(source))
      if (present(references)) ierr = nf90_put_att(ncid, NF90_GLOBAL, 'references',  trim(references))
      if (present(comment))    ierr = nf90_put_att(ncid, NF90_GLOBAL, 'comment',     trim(comment))
      if (present(ncdata))     ierr = nf90_put_att(ncid, NF90_GLOBAL, 'ncdata',      trim(ncdata))
      if (present(calculation_type)) ierr = nf90_put_att(ncid, NF90_GLOBAL, &
                                            'calculation_type',  trim(calculation_type))

      ! CF hint (pick the version you follow; or omit if you don’t want to claim one)
      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.8')

      call check_nc(ierr,'put_att(globals)')
      ierr = nf90_enddef(ncid); call check_nc(ierr,'enddef(globals)')
    end subroutine ncfile_set_globals

   !----------------------------
   ! write a physically 1D field with dims ({level,ilevel})
   subroutine ncfile_put_col1d_notime(varname, field, units, long_name)
      character(len=*), intent(in) :: varname
      real(r8),         intent(in) :: field(:)
      character(len=*), intent(in), optional :: units, long_name

      integer :: ierr, varid, dimid_zX, nzX_save
      integer :: dimids(1)
      integer :: start(1), count(1)


      ! Determine vertical dim (lev or ilev) and check sizes 
      if (size(field,1) == nz_save) then
         dimid_zX = dimid_z
         nzX_save = nz_save
      elseif (size(field,1) == nze_save) then
         dimid_zX = dimid_ze
         nzX_save = nze_save
      else
         write(*,*) 'Size mismatch in ncfile_put_col3d for ', trim(varname)
         stop 1
      end if

      ! Try to find variable
      ierr = nf90_inq_varid(ncid, trim(varname), varid)

      if (ierr /= nf90_noerr) then
         ! Need to (re)enter define mode and create it
         ierr = nf90_redef(ncid); call check_nc(ierr, 'nf90_redef')

         dimids = (/ dimid_zX /)
         ierr   = nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varid)
         call check_nc(ierr, 'def_var '//trim(varname))

         if (present(units)) then
            ierr = nf90_put_att(ncid, varid, 'units', units)
            call check_nc(ierr, 'put_att units '//trim(varname))
         end if

         if (present(long_name)) then
            ierr = nf90_put_att(ncid, varid, 'long_name', long_name)
            call check_nc(ierr, 'put_att long_name '//trim(varname))
         end if

         ierr = nf90_enddef(ncid)
         call check_nc(ierr, 'nf90_enddef(new var)')
      end if

      ! Now write this time slice
      start = (/ 1 /)
      count = (/ nzX_save /)

      ierr = nf90_put_var(ncid, varid, field, start=start, count=count)
      call check_nc(ierr, 'put_var '//trim(varname))


      
    end subroutine ncfile_put_col1d_notime

   !----------------------------
   ! write a 3D field with dims (lon,lat,level,time)
   subroutine ncfile_put_3d(varname, field, itime, units, long_name)
      character(len=*), intent(in) :: varname
      real,             intent(in) :: field(:,:,:)
      integer,          intent(in) :: itime  ! 1-based time index
      character(len=*), intent(in), optional :: units, long_name

      integer :: ierr, varid
      integer :: dimids(4)
      integer :: start(4), count(4)

      ! Sanity check on sizes (optional but helpful)
      if (size(field,1) /= ncol_save .or. &
          size(field,2) /= nz_save) then
         write(*,*) 'Size mismatch in ncfile_put_3d for ', trim(varname)
         stop 1
      end if

      ! Try to find variable
      ierr = nf90_inq_varid(ncid, trim(varname), varid)

      if (ierr /= nf90_noerr) then
         ! Need to (re)enter define mode and create it
         ierr = nf90_redef(ncid); call check_nc(ierr, 'nf90_redef')

         dimids = (/ dimid_x, dimid_y, dimid_z, dimid_time /)
         ierr   = nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varid)
         call check_nc(ierr, 'def_var '//trim(varname))

         if (present(units)) then
            ierr = nf90_put_att(ncid, varid, 'units', units)
            call check_nc(ierr, 'put_att units '//trim(varname))
         end if

         if (present(long_name)) then
            ierr = nf90_put_att(ncid, varid, 'long_name', long_name)
            call check_nc(ierr, 'put_att long_name '//trim(varname))
         end if

         ierr = nf90_enddef(ncid)
         call check_nc(ierr, 'nf90_enddef(new var)')
      end if

      ! Now write this time slice
      start = (/ 1, 1, 1, itime /)
      count = (/ nx_save, ny_save, nz_save, 1 /)

      ierr = nf90_put_var(ncid, varid, field, start=start, count=count)
      call check_nc(ierr, 'put_var '//trim(varname))
   end subroutine ncfile_put_3d

   !----------------------------
   ! write a physically 2D field with dims (col), i.e., no time dim ...
   subroutine ncfile_put_col2d_notime(varname, field, units, long_name)
      character(len=*), intent(in) :: varname
      real(r8),         intent(in) :: field(:)
      character(len=*), intent(in), optional :: units, long_name

      integer :: ierr, varid
      integer :: dimids(1)
      integer :: start(1), count(1)

      ! Sanity check on sizes (optional but helpful)
      if (size(field,1) /= ncol_save ) then
         write(*,*) 'Size mismatch in ncfile_put_col2d for ', trim(varname)
         stop 1
      end if

      ! Try to find variable
      ierr = nf90_inq_varid(ncid, trim(varname), varid)

      if (ierr /= nf90_noerr) then
         ! Need to (re)enter define mode and create it
         ierr = nf90_redef(ncid); call check_nc(ierr, 'nf90_redef')

         dimids = (/ dimid_col /)
         ierr   = nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varid)
         call check_nc(ierr, 'def_var '//trim(varname))

         if (present(units)) then
            ierr = nf90_put_att(ncid, varid, 'units', units)
            call check_nc(ierr, 'put_att units '//trim(varname))
         end if

         if (present(long_name)) then
            ierr = nf90_put_att(ncid, varid, 'long_name', long_name)
            call check_nc(ierr, 'put_att long_name '//trim(varname))
         end if

         ierr = nf90_enddef(ncid)
         call check_nc(ierr, 'nf90_enddef(new var)')
      end if

      ! Now write this time slice
      start = (/ 1 /)
      count = (/ ncol_save /)

      ierr = nf90_put_var(ncid, varid, field, start=start, count=count)
      call check_nc(ierr, 'put_var '//trim(varname))
    end subroutine ncfile_put_col2d_notime

   !----------------------------
   ! write a physically 2D field with dims (col,time)
   subroutine ncfile_put_col2d(varname, field, itime, units, long_name)
      character(len=*), intent(in) :: varname
      real(r8),         intent(in) :: field(:)
      integer,          intent(in) :: itime  ! 1-based time index
      character(len=*), intent(in), optional :: units, long_name

      integer :: ierr, varid
      integer :: dimids(2)
      integer :: start(2), count(2)

      ! Sanity check on sizes (optional but helpful)
      if (size(field,1) /= ncol_save ) then
         write(*,*) 'Size mismatch in ncfile_put_col2d for ', trim(varname)
         stop 1
      end if

      ! Try to find variable
      ierr = nf90_inq_varid(ncid, trim(varname), varid)

      if (ierr /= nf90_noerr) then
         ! Need to (re)enter define mode and create it
         ierr = nf90_redef(ncid); call check_nc(ierr, 'nf90_redef')

         dimids = (/ dimid_col, dimid_time /)
         ierr   = nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varid)
         call check_nc(ierr, 'def_var '//trim(varname))

         if (present(units)) then
            ierr = nf90_put_att(ncid, varid, 'units', units)
            call check_nc(ierr, 'put_att units '//trim(varname))
         end if

         if (present(long_name)) then
            ierr = nf90_put_att(ncid, varid, 'long_name', long_name)
            call check_nc(ierr, 'put_att long_name '//trim(varname))
         end if

         ierr = nf90_enddef(ncid)
         call check_nc(ierr, 'nf90_enddef(new var)')
      end if

      ! Now write this time slice
      start = (/ 1, itime /)
      count = (/ ncol_save, 1 /)

      ierr = nf90_put_var(ncid, varid, field, start=start, count=count)
      call check_nc(ierr, 'put_var '//trim(varname))
    end subroutine ncfile_put_col2d

   !----------------------------
   ! write a physically 3D field with dims (col,level,time)
   subroutine ncfile_put_col3d(varname, field, itime, units, long_name)
      character(len=*), intent(in) :: varname
      real(r8),         intent(in) :: field(:,:)
      integer,          intent(in) :: itime  ! 1-based time index
      character(len=*), intent(in), optional :: units, long_name

      integer :: ierr, varid, dimid_zX, nzX_save
      integer :: dimids(3)
      integer :: start(3), count(3)

      ! Determine vertical dim (lev or ilev) and check sizes 
      if (size(field,2) == nz_save) then
         dimid_zX = dimid_z
         nzX_save = nz_save
      elseif (size(field,2) == nze_save) then
         dimid_zX = dimid_ze
         nzX_save = nze_save
      else
         write(*,*) 'Size mismatch in ncfile_put_col3d for ', trim(varname)
         stop 1
      end if
      if (size(field,1) /= ncol_save) then
         write(*,*) 'Size mismatch in ncfile_put_col3d for ', trim(varname)
         stop 1
      end if

      ! Try to find variable
      ierr = nf90_inq_varid(ncid, trim(varname), varid)

      if (ierr /= nf90_noerr) then
         ! Need to (re)enter define mode and create it
         ierr = nf90_redef(ncid); call check_nc(ierr, 'nf90_redef')

         dimids = (/ dimid_col, dimid_zX , dimid_time /)
         ierr   = nf90_def_var(ncid, trim(varname), NF90_REAL, dimids, varid)
         call check_nc(ierr, 'def_var '//trim(varname))

         if (present(units)) then
            ierr = nf90_put_att(ncid, varid, 'units', units)
            call check_nc(ierr, 'put_att units '//trim(varname))
         end if

         if (present(long_name)) then
            ierr = nf90_put_att(ncid, varid, 'long_name', long_name)
            call check_nc(ierr, 'put_att long_name '//trim(varname))
         end if

         ierr = nf90_enddef(ncid)
         call check_nc(ierr, 'nf90_enddef(new var)')
      end if

      ! Now write this time slice
      start = (/ 1, 1, itime /)
      count = (/ ncol_save, nzX_save, 1 /)

      ierr = nf90_put_var(ncid, varid, field, start=start, count=count)
      call check_nc(ierr, 'put_var '//trim(varname))

    end subroutine ncfile_put_col3d

   !----------------------------
   subroutine ncfile_close()
      integer :: ierr
      if (ncid /= -1) then
         ierr = nf90_close(ncid)
         call check_nc(ierr, 'nf90_close')
         ncid = -1
      end if
   end subroutine ncfile_close

end module nc_flexout_mod
