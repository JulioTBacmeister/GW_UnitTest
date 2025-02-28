! !MODULE: shr_string_mod -- string and list methods
!
! !DESCRIPTION:
!    General string and specific list method.  A list is a single string
!    that is delimited by a character forming multiple fields, ie,
!    character(len=*) :: mylist = "t:s:u1:v1:u2:v2:taux:tauy"
!    The delimiter is called listDel in this module, is default ":",
!    but can be set by a call to shr_string_listSetDel.
!
! !REVISION HISTORY:
!     2005-Apr-28 - T. Craig - first version
!     2025-Feb-22 - J. Bacmeister ;) totally hacked to work in GW unit test
!
! !INTERFACE: ------------------------------------------------------------------

module shr_string_mod

  ! !USES:
!++jtb (better not need this)
!!! #include "shr_assert.h"
  use shr_kind_mod   ! F90 kinds

  
  !++jtb
  ! use shr_sys_mod    ! shared system calls
  ! use shr_timer_mod, only : shr_timer_get, shr_timer_start, shr_timer_stop
  ! use shr_log_mod,   only : errMsg    => shr_log_errMsg
  ! use shr_log_mod,   only : s_loglev  => shr_log_Level
  ! use shr_log_mod,   only : s_logunit => shr_log_Unit

  implicit none
  private

  ! !PUBLIC TYPES:

  ! no public types

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_string_countChar       ! Count number of char in string, fn
  public :: shr_string_toUpper         ! Convert string to upper-case
  public :: shr_string_toLower         ! Convert string to lower-case

  ! !PUBLIC DATA MEMBERS:

  ! no public data members

  !EOP

  character(len=1)    ,save :: listDel  = ":"    ! note single exec implications
  character(len=2)    ,save :: listDel2 = "::"   ! note single exec implications
  logical             ,save :: doabort  = .true.
  integer(SHR_KIND_IN),save :: debug    = 0

  !===============================================================================
contains
  !===============================================================================

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_string_countChar -- Count number of occurances of a character
  !
  ! !DESCRIPTION:
  !  count number of occurances of a single character in a string
  !     \newline
  !     n = shr\_string\_countChar(string,character)
  !
  ! !REVISION HISTORY:
  !     2005-Feb-28 - First version from dshr_bundle
  !
  ! !INTERFACE: ------------------------------------------------------------------

  integer function shr_string_countChar(str,char,rc)


    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)        ,intent(in)           :: str   ! string to search
    character(1)        ,intent(in)           :: char  ! char to search for
    integer(SHR_KIND_IN),intent(out),optional :: rc    ! return code

    !EOP

    !----- local -----
    integer(SHR_KIND_IN) :: count    ! counts occurances of char
    integer(SHR_KIND_IN) :: n        ! generic index
    integer(SHR_KIND_IN) :: t01 = 0  ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_countChar) "
    character(*),parameter :: F00     = "('(shr_string_countChar) ',4a)"

    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    !++jtb
    !if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    !if (debug>1) call shr_timer_start(t01)

    count = 0
    do n = 1, len_trim(str)
       if (str(n:n) == char) count = count + 1
    end do
    shr_string_countChar = count

    if (present(rc)) rc = 0

    !++jtb
    !if (debug>1) call shr_timer_stop (t01)

  end function shr_string_countChar

  !===============================================================================
  !BOP ===========================================================================
  ! !IROUTINE: shr_string_toUpper -- Convert string to upper case
  !
  ! !DESCRIPTION:
  !     Convert the input string to upper-case.
  !     Use achar and iachar intrinsics to ensure use of ascii collating sequence.
  !
  ! !REVISION HISTORY:
  !     2005-Dec-20 - Move CAM version over to shared code.
  !
  ! !INTERFACE: ------------------------------------------------------------------

  function shr_string_toUpper(str)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: str      ! String to convert to upper case
    character(len=len(str))      :: shr_string_toUpper

    !----- local -----
    integer(SHR_KIND_IN) :: i             ! Index
    integer(SHR_KIND_IN) :: aseq          ! ascii collating sequence
    integer(SHR_KIND_IN) :: LowerToUpper  ! integer to convert case
    character(len=1)     :: ctmp          ! Character temporary
    integer(SHR_KIND_IN) :: t01 = 0       ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_toUpper) "
    character(*),parameter :: F00     = "('(shr_string_toUpper) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !++jtb
    !if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    !if (debug>1) call shr_timer_start(t01)

    LowerToUpper = iachar("A") - iachar("a")

    do i = 1, len(str)
       ctmp = str(i:i)
       aseq = iachar(ctmp)
       if ( aseq >= iachar("a") .and. aseq <= iachar("z") ) &
            ctmp = achar(aseq + LowertoUpper)
       shr_string_toUpper(i:i) = ctmp
    end do

    !++jtb
    !if (debug>1) call shr_timer_stop (t01)

  end function shr_string_toUpper

  !===============================================================================
  !BOP ===========================================================================
  ! !IROUTINE: shr_string_toLower -- Convert string to lower case
  !
  ! !DESCRIPTION:
  !     Convert the input string to lower-case.
  !     Use achar and iachar intrinsics to ensure use of ascii collating sequence.
  !
  ! !REVISION HISTORY:
  !     2006-Apr-20 - Creation
  !
  ! !INTERFACE: ------------------------------------------------------------------
  function shr_string_toLower(str)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: str      ! String to convert to lower case
    character(len=len(str))      :: shr_string_toLower

    !----- local -----
    integer(SHR_KIND_IN) :: i            ! Index
    integer(SHR_KIND_IN) :: aseq         ! ascii collating sequence
    integer(SHR_KIND_IN) :: UpperToLower ! integer to convert case
    character(len=1)     :: ctmp         ! Character temporary
    integer(SHR_KIND_IN) :: t01 = 0      ! timer

    !----- formats -----
    character(*),parameter :: subName =   "(shr_string_toLower) "
    character(*),parameter :: F00     = "('(shr_string_toLower) ',4a)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !++jtb 
    !if (debug>1 .and. t01<1) call shr_timer_get(t01,subName)
    !if (debug>1) call shr_timer_start(t01)

    UpperToLower = iachar("a") - iachar("A")

    do i = 1, len(str)
       ctmp = str(i:i)
       aseq = iachar(ctmp)
       if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
            ctmp = achar(aseq + UpperToLower)
       shr_string_toLower(i:i) = ctmp
    end do

    !++jtb
    !if (debug>1) call shr_timer_stop (t01)

  end function shr_string_toLower

  !===============================================================================
  !===============================================================================

end module shr_string_mod
