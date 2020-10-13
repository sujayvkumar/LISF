!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LISroutingWL_obsMod
!BOP
! 
! !MODULE: LISroutingWL_obsMod
! 
! !DESCRIPTION: 
!  This module handles the use of a LIS model simulation 
!  output as "observations" for data assimilation. This
!  plugin is typically used to handle the computations of 
!  scaling factors such as cumulative distribution function
!  (CDF) for use in DA
! 
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
!

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LISroutingWL_obsInit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: routingwlobs
!
!EOP
  
  type, public :: routingwlobsdec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr
     real          :: datares
     real          :: run_dd(8)
     character*50  :: map_proj
     character*50  :: format
     character*50  :: wstyle
     character*50  :: wopt
     character*100 :: odir
     character*20  :: security_class
     character*20  :: distribution_class
     character*20  :: data_category
     character*20  :: area_of_data
     character*20  :: write_interval
!--------------------------------------------------------
!  interpolation/upscaling weights
!--------------------------------------------------------
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,  allocatable     :: w11(:)
     real,  allocatable     :: w12(:)
     real,  allocatable     :: w21(:)
     real,  allocatable     :: w22(:)

  end type routingwlobsdec

  type(routingwlobsdec)  :: routingwlobs

contains

!BOP
! !ROUTINE: LISroutingWL_obsInit
! \label{LISroutingWL_obsInit}
! 
! !INTERFACE: 
  subroutine LISroutingWL_obsInit()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_logMod

    implicit none
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of a 
! land surface model output (from a LIS simulation) as observations.  
! 
!EOP
    integer                 :: n
    integer                 :: rc
    real                    :: gridDesci(20)

    n = 1

    routingwlobs%run_dd             = LDT_rc%udef

    routingwlobs%security_class     = ''
    routingwlobs%distribution_class = ''
    routingwlobs%data_category      = ''
    routingwlobs%area_of_data       = ''
    routingwlobs%write_interval     = ''

    call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%format, &
         label="LIS water level output format:",rc=rc)
    call LDT_verify(rc,'LIS water level output format: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%wopt, &
         label="LIS water level output methodology:",rc=rc)
    call LDT_verify(rc,'LIS water level output methodology: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%wstyle, &
         label="LIS water level output naming style:",rc=rc)
    call LDT_verify(rc,'LIS water level output naming style: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%map_proj, &
         label="LIS water level output map projection:",rc=rc)
    call LDT_verify(rc,'LIS water level output map projection: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%nest, &
         label="LIS water level output nest index:",rc=rc)
    call LDT_verify(rc,'LIS water level output nest index: not defined')

    call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%odir, &
         label="LIS water level output directory:",rc=rc)
    call LDT_verify(rc,'LIS water level output directory: not defined')

    ! WMO-convention specific identifiers
    if ( routingwlobs%wstyle == "WMO convention") then 
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%security_class, &
       label="LIS water level security class:",rc=rc)
       call LDT_verify(rc,'LIS water level security class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%distribution_class, &
       label="LIS water level distribution class:",rc=rc)
       call LDT_verify(rc,'LIS water level distribution class: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%data_category, &
       label="LIS water level data category:",rc=rc)
       call LDT_verify(rc,'LIS water level data category: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%area_of_data, &
       label="LIS water level area of data:",rc=rc)
       call LDT_verify(rc,'LIS water level area of data: not defined')

       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%write_interval, &
       label="LIS water level write interval:",rc=rc)
       call LDT_verify(rc,'LIS water level write interval: not defined')
    endif

    if(routingwlobs%map_proj.eq."latlon") then 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(1),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(2),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain upper right lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(3),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain upper right lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(4),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain resolution (dx):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(5),rc=rc)
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain resolution (dy):",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(6),rc=rc)       
       
       routingwlobs%datares = min(routingwlobs%run_dd(5),routingwlobs%run_dd(6))

       routingwlobs%nc    = (nint((routingwlobs%run_dd(4)-routingwlobs%run_dd(2))/&
            routingwlobs%run_dd(5))) + 1
       routingwlobs%nr    = (nint((routingwlobs%run_dd(3)-routingwlobs%run_dd(1))/&
            routingwlobs%run_dd(6))) + 1

       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = routingwlobs%nc
       gridDesci(3) = routingwlobs%nr
       gridDesci(4) = routingwlobs%run_dd(1)
       gridDesci(5) = routingwlobs%run_dd(2)
       gridDesci(6) = 128
       gridDesci(7) = routingwlobs%run_dd(3)
       gridDesci(8) = routingwlobs%run_dd(4)
       gridDesci(9) = routingwlobs%run_dd(5)
       gridDesci(10) = routingwlobs%run_dd(6)
       gridDesci(20) = 64

    elseif(routingwlobs%map_proj.eq."lambert") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(1),rc=rc)
       call LDT_verify(rc,'LIS water level domain lower left lat: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(2),rc=rc)
       call LDT_verify(rc,'LIS water level domain lower left lon: not defined')
       
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(3),rc=rc)
       call LDT_verify(rc,'LIS water level domain true lat1: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain true lat2:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(4),rc=rc)
       call LDT_verify(rc,'LIS water level domain true lat2: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(5),rc=rc)
       call LDT_verify(rc,'LIS water level domain standard lon: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(6),rc=rc)
       call LDT_verify(rc,'LIS water level domain resolution: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(7),rc=rc)
       call LDT_verify(rc,'LIS water level domain x-dimension size: not defined')
       
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(8),rc=rc)
       call LDT_verify(rc,'LIS water level domain y-dimension size: not defined')

       routingwlobs%datares = routingwlobs%run_dd(6)/100.0

       routingwlobs%nc    = routingwlobs%run_dd(7)
       routingwlobs%nr    = routingwlobs%run_dd(8)

       gridDesci = 0
       gridDesci(1) = 3 
       gridDesci(2) = routingwlobs%nc
       gridDesci(3) = routingwlobs%nr
       gridDesci(4) = routingwlobs%run_dd(1)
       gridDesci(5) = routingwlobs%run_dd(2)
       gridDesci(6) = 8
       gridDesci(7) = routingwlobs%run_dd(4)
       gridDesci(8) = routingwlobs%run_dd(6)
       gridDesci(9) = routingwlobs%run_dd(6)
       gridDesci(10) = routingwlobs%run_dd(3)
       gridDesci(11) = routingwlobs%run_dd(5)
       gridDesci(20) = 64

    elseif(routingwlobs%map_proj.eq."polar") then 
       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain lower left lat:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(1),rc=rc)
 

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain lower left lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(2),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain true lat1:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(3),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain true lat2:",rc=rc)       
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(4),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain standard lon:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(5),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain resolution:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(6),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain x-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(7),rc=rc)

       call ESMF_ConfigFindLabel(LDT_config,&
            "LIS water level domain y-dimension size:",rc=rc)
       call ESMF_ConfigGetAttribute(LDT_config,routingwlobs%run_dd(8),rc=rc)

    endif

!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------
    if(LDT_isLDTatAfinerResolution(n,routingwlobs%datares)) then 

       allocate(routingwlobs%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(routingwlobs%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(routingwlobs%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(routingwlobs%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
      
       allocate(routingwlobs%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(routingwlobs%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(routingwlobs%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       allocate(routingwlobs%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call bilinear_interp_input(n, gridDesci, &
            routingwlobs%n11, &
            routingwlobs%n12, routingwlobs%n21, &
            routingwlobs%n22, routingwlobs%w11, &
            routingwlobs%w12, routingwlobs%w21, &
            routingwlobs%w22)

    else

       allocate(routingwlobs%n11(&
            routingwlobs%nc*&
            routingwlobs%nr))

       call upscaleByAveraging_input(&
            gridDesci,&
            LDT_rc%gridDesc(n,:),&
            routingwlobs%nc*routingwlobs%nr,&
            LDT_rc%lnc(n)*LDT_rc%lnr(n),&
            routingwlobs%n11)
       
    endif
    
!  which variable we want in the DA obs computations. 
    call LDT_initializeDAobsEntry(LDT_DAobsData(1)%wl_obs, "m3/m3",1,1)
    LDT_DAobsData(1)%wl_obs%selectStats = 1

  end subroutine LISroutingWL_obsInit
  
end module LISroutingWL_obsMod
