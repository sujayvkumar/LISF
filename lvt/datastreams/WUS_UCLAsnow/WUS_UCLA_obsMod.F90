!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: WUS_UCLA_obsMod
! \label(WUS_UCLA_obsMod)
!
! !INTERFACE:
module WUS_UCLA_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 May 2010   Sujay Kumar  Initial Specification
! 
!EOP
! !USES: 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: WUS_UCLA_obsinit !Initializes structures for reading WUS_UCLA data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: WUSUCLAobs !Object to hold WUS_UCLA observation attributes
!EOP
  type, public :: wusuclaobsdec
     character*100           :: odir
     integer                 :: nc
     integer                 :: nr
     integer                 :: c_off, r_off
     real                    :: udef
     real,  allocatable      :: snwd(:,:)
     integer, allocatable    :: n11(:)
     logical                 :: startFlag

  end type wusuclaobsdec

  type(wusuclaobsdec), allocatable :: wusuclaobs(:)

contains
  
!BOP
! 
! !ROUTINE: WUS_UCLA_obsInit
! \label{WUS_UCLA_obsInit}
!
! !INTERFACE:
  subroutine WUS_UCLA_obsinit(i)
! 
! !USES:   
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
!  for reading WUS_UCLA data. The WUS_UCLA data is provides in the WGS 1984 grid
!  with 30 arc second resolution. The domain extents are from (24.9504, -124.7337) 
!  to (52.8754, -66.9421). The data entries are 16-bit signed integers. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: status, rc
    integer            :: ftn, k
    real*8             :: tdur
    integer            :: syr, smo, sda, shr, smn, sss
    integer            :: eyr, emo, eda, ehr, emn, ess
    integer            :: ts
    character*100      :: coordfile
    character*100      :: mdata
    real               :: xi1,xj1,xmesh,orient,alat1,alon1
    integer            :: t
    real               :: gridDesci(50)
    real                   :: cornerlat1, cornerlat2
    real                   :: cornerlon1, cornerlon2    

    
    if(.not.allocated(wusuclaobs)) then 
       allocate(wusuclaobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, wusuclaobs(i)%odir, &
         label='WUS UCLA observation directory:',rc=status)
    call LVT_verify(status, 'WUS_UCLA observation directory: not defined')

    wusuclaobs(i)%udef = -9999.0
    wusuclaobs(i)%startFlag = .true. 
    
    ts = 86400
    call LVT_update_timestep(LVT_rc, 86400)


    cornerlat1 = max(31.002222, &
         nint((LVT_rc%gridDesc(4)-31.002222)/&
         0.0044444)*0.0044444+31.002222-50*0.0044444)
    cornerlon1 = max(-124.99777777778, &
         nint((LVT_rc%gridDesc(5)+&
         124.997777778)/0.0044444)*0.0044444-124.99777777778-50*0.0044444)
    cornerlat2 = min(48.997777778, &
         nint((LVT_rc%gridDesc(7)-31.002222)/&
         0.0044444)*0.0044444+31.002222+50*0.0044444)
    cornerlon2 = min(-102.002222222222, &
         nint((LVT_rc%gridDesc(8)+124.997777778)/&
         0.0044444)*0.0044444-124.99777777778+50*0.0044444)

    wusuclaobs(i)%c_off = nint((cornerlon1 + 124.997778)/0.0044444)+1
    wusuclaobs(i)%r_off = nint((cornerlat1 - 31.002222)/0.0044444)+1
    
    
    wusuclaobs(i)%nc = nint((cornerlon2-cornerlon1)/0.0044444)+1
    wusuclaobs(i)%nr = nint((cornerlat2-cornerlat1)/0.0044444)+1
       

    gridDesci = 0 
    gridDesci(1) = 0
    gridDesci(2) = wusuclaobs(i)%nc
    gridDesci(3) = wusuclaobs(i)%nr
    griddesci(4) = cornerlat1
    gridDesci(5) = cornerlon1
    gridDesci(6) = 128
    gridDesci(7) = cornerlat2
    gridDesci(8) = cornerlon2
    gridDesci(9) = 0.0044444
    gridDesci(10) =0.0044444
    gridDesci(20) = 64


    allocate(wusuclaobs(i)%n11(wusuclaobs(i)%nc*wusuclaobs(i)%nr))


    call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc,&
         wusuclaobs(i)%nc*wusuclaobs(i)%nr, &
         LVT_rc%lnc*LVT_rc%lnr, wusuclaobs(i)%n11)
    

    allocate(wusuclaobs(i)%snwd(LVT_rc%lnc,LVT_rc%lnr))
    wusuclaobs(i)%snwd = LVT_rc%udef

  end subroutine WUS_UCLA_obsinit


end module WUS_UCLA_obsMod
