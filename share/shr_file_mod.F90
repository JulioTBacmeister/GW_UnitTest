! !MODULE: shr_file_mod.F90 --- Module to handle various file utilily functions.
!
! !DESCRIPTION:
!
! Miscilaneous methods to handle file and directory utilities as well as FORTRAN
! unit control. Also put/get local files into/from archival location
!
! File utilites used with CCSM Message passing:
!
! shr_file_stdio is the main example here, it changes the working directory,
!                changes stdin and stdout to a given filename.
!
! This is needed because some implementations of MPI with MPMD so that
! each executable can run in a different working directory and redirect
! output to different files.
!
! File name archival convention, eg.
!    call shr_file_put(rcode,"foo","mss:/USER/foo",rtpd=3650)
! is extensible -- the existence of the option file name prefix, eg. "mss:",
! and optional arguments, eg. rtpd-3650 can be used to access site-specific
! storage devices.  Based on CCM (atmosphere) getfile & putfile routines, but
! intended to be a more extensible, shared code.
!
! !REVISION HISTORY:
!   2006-05-08 E. Kluzek, Add in shr_file_mod and getUnit, freeUnif methods.
!   2000-??-?? B. Kauffman, original version circa 2000
!
! !INTERFACE: ------------------------------------------------------------------

MODULE shr_file_mod

  ! !USES:

  use shr_kind_mod  ! defines kinds
  !++jtb
  use, intrinsic :: iso_fortran_env, only: output_unit

  IMPLICIT none

  PRIVATE           ! By default everything is private to this module

  ! !PUBLIC TYPES:

  ! no public types

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_file_getUnit      ! Get a logical unit for reading or writing
  public :: shr_file_freeUnit     ! Free a logical unit

  ! !PUBLIC DATA MEMBERS:

  ! Integer flags for recognized prefixes on file get/put operations
  integer(SHR_KIND_IN), parameter, public :: shr_file_noPrefix      = 0 ! no recognized prefix
  integer(SHR_KIND_IN), parameter, public :: shr_file_nullPrefix    = 1 ! null:
  integer(SHR_KIND_IN), parameter, public :: shr_file_cpPrefix      = 2 ! cp:
  integer(SHR_KIND_IN), parameter, public :: shr_file_mssPrefix     = 3 ! mss:
  integer(SHR_KIND_IN), parameter, public :: shr_file_hpssPrefix    = 4 ! hpss:

  !++jtb
  ! low-level shared variables for logging, these may not be parameters
  integer(SHR_KIND_IN) :: s_loglev = 0
  integer(SHR_KIND_IN) :: s_logunit  = output_unit

  !EOP
  !--- unit numbers, users can ask for unit numbers from 0 to min, but getUnit
  !--- won't give a unit below min, users cannot ask for unit number above max
  !--- for backward compatability.
  !--- eventually, recommend min as hard lower limit (tcraig, 9/2007)
  integer(SHR_KIND_IN),parameter :: shr_file_minUnit = 10      ! Min unit number to give
  integer(SHR_KIND_IN),parameter :: shr_file_maxUnit = 99      ! Max unit number to give
  logical, save :: UnitTag(0:shr_file_maxUnit) = .false. ! Logical units in use

  !===============================================================================
CONTAINS
  !===============================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_file_getUnit -- Get a free FORTRAN unit number
  !
  ! !DESCRIPTION: Get the next free FORTRAN unit number.
  !
  ! !REVISION HISTORY:
  !     2005-Dec-14 - E. Kluzek - creation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  INTEGER FUNCTION shr_file_getUnit ( unit )

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in),optional :: unit ! desired unit number

    !EOP

    !----- local -----
    integer(SHR_KIND_IN)   :: n      ! loop index
    logical                :: opened ! If unit opened or not

    !----- formats -----
    character(*),parameter :: subName = '(shr_file_getUnit) '
    character(*),parameter :: F00   = "('(shr_file_getUnit) ',A,I4,A)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------
    shr_file_getUnit = -1
    if (present (unit)) then
       inquire( unit, opened=opened )
       if (unit < 0 .or. unit > shr_file_maxUnit) then
          write(s_logunit,F00) 'invalid unit number request:', unit
          !++jtb (disable this replace w/ stop although really wgaF )
          !call shr_sys_abort( 'ERROR: bad input unit number' )
          STOP
       else if (opened .or. UnitTag(unit) .or. unit == 0 .or. unit == 5 &
            .or. unit == 6) then
          write(s_logunit,F00) 'unit number ', unit, ' is already in use'
          !++jtb (disable this replace w/ stop although really wgaF )
          !call shr_sys_abort( 'ERROR: Input unit number already in use' )
          STOP
       else
          shr_file_getUnit = unit
          UnitTag (unit)   = .true.
          return
       end if

    else
       ! --- Choose first available unit other than 0, 5, or 6  ------
       do n=shr_file_maxUnit, shr_file_minUnit, -1
          inquire( n, opened=opened )
          if (n == 5 .or. n == 6 .or. opened) then
             cycle
          end if
          if ( .not. UnitTag(n) ) then
             shr_file_getUnit = n
             UnitTag(n)       = .true.
             return
          end if
       end do
    end if

    !++jtb (disable this replace w/ stop although really wgaF )
    !call shr_sys_abort( subName//': Error: no available units found' )
    STOP
    
  END FUNCTION shr_file_getUnit

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_file_freeUnit -- Free up a FORTRAN unit number
  !
  ! !DESCRIPTION: Free up the given unit number
  !
  ! !REVISION HISTORY:
  !     2005-Dec-14 - E. Kluzek - creation
  !
  ! !INTERFACE: ------------------------------------------------------------------

  SUBROUTINE shr_file_freeUnit ( unit)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in) :: unit  ! unit number to be freed

    !EOP

    !----- local -----

    !----- formats -----
    character(*), parameter :: subName = '(shr_file_freeUnit) '
    character(*), parameter :: F00 =   "('(shr_file_freeUnit) ',A,I4,A)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (unit < 0 .or. unit > shr_file_maxUnit) then
       if (s_loglev > 0) write(s_logunit,F00) 'invalid unit number request:', unit
    else if (unit == 0 .or. unit == 5 .or. unit == 6) then
       !++jtb (disable this replace w/ stop although really wgaF )
       !call shr_sys_abort( subName//': Error: units 0, 5, and 6 must not be freed' )
       STOP
    else if (UnitTag(unit)) then
       UnitTag (unit) = .false.
    else
       if (s_loglev > 0) write(s_logunit,F00) 'unit ', unit, ' was not in use'
    end if

    return

  END SUBROUTINE shr_file_freeUnit


  !===============================================================================

END MODULE shr_file_mod
