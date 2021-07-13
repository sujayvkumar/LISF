!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: NOAAGLDAS_obsMod
! \label(NOAAGLDAS_obsMod)
!
! !INTERFACE:
module NOAAGLDAS_obsMod
! 
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the
!  NOAA GLDAS output
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 July 2021: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NOAAGLDAS_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NOAAGLDASobs
!EOP
  type, public :: lvtnoaaobsdec
     character*100              :: odir
     real                       :: gridDesci(50)
     integer                    :: nc, nr
     type(proj_info)            :: proj
     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)
     real,  allocatable         :: w11(:)
     real,  allocatable         :: w12(:)
     real,  allocatable         :: w21(:)
     real,  allocatable         :: w22(:)
     real,  allocatable         :: rlat(:)
     real,  allocatable         :: rlon(:)
     real,    allocatable       :: smobs(:,:)
     logical                    :: startmode     
     type(ESMF_Time)            :: startTime
     type(ESMF_TimeInterval)    :: dt
  end type lvtnoaaobsdec

  type(lvtnoaaobsdec), allocatable:: NOAAGLDASobs(:)

contains
  
!BOP
! 
! !ROUTINE: NOAAGLDAS_obsInit
! \label{NOAAGLDAS_obsInit}
!
! !INTERFACE: 
  subroutine NOAAGLDAS_obsinit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading ESACCI AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    if(.not.allocated(NOAAGLDASobs)) then 
       allocate(NOAAGLDASobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, NOAAGLDASobs(i)%odir, &
         label='NOAA GLDAS data directory:', rc=status)
    call LVT_verify(status, 'NOAA GLDAS data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    NOAAGLDASobs(i)%startmode = .true. 

    allocate(NOAAGLDASobs(i)%smobs(LVT_rc%lnc,LVT_rc%lnr))

    NOAAGLDASobs(i)%nc = 1536
    NOAAGLDASobs(i)%nr = 768
    
    gridDesci = 0 
    gridDesci(1) = 4 
    gridDesci(2) = 1536
    gridDesci(3) = 768
    gridDesci(4) = 89.8207074233875
    gridDesci(5) = 0 
    gridDesci(6) = 128
    gridDesci(7) = -89.8207074233875
    gridDesci(8) = -0.234375 ! 359.765625
    gridDesci(9) =  0.234375
    gridDesci(10) = 384
    gridDesci(20) = 64
    
    allocate(NOAAGLDASobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(NOAAGLDASobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(NOAAGLDASobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(NOAAGLDASobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,NOAAGLDASobs(i)%rlat, &
         NOAAGLDASobs(i)%rlon,NOAAGLDASobs(i)%n11, &
         NOAAGLDASobs(i)%n12, NOAAGLDASobs(i)%n21, &
         NOAAGLDASobs(i)%n22, NOAAGLDASobs(i)%w11, &
         NOAAGLDASobs(i)%w12, NOAAGLDASobs(i)%w21, &
         NOAAGLDASobs(i)%w22)


    call ESMF_TimeSet(NOAAGLDASobs(i)%startTime,&
         yy=2012,&
         mm = 7,&
         dd = 1,&
         h = 0, &
         m = 0, &
         s = 0, &
         calendar=LVT_calendar, &
         rc=status)
    call ESMF_TimeIntervalSet(NOAAGLDASobs(i)%dt, s=86400,rc=status)
    
  end subroutine NOAAGLDAS_obsinit


end module NOAAGLDAS_obsMod
