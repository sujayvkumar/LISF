!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: LPRMvod_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! LPRM vegetation optical depth (vod) retrievals
!
! !REVISION HISTORY: 
!  28 May 2019: Sujay Kumar, Initial Specification
!
#include "LDT_misc.h"
module hydrowebWL_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: hydrowebWL_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: hydrowebWLobs
!EOP
  type, public :: hydrowebwlobsdec
     type(ESMF_TimeInterval) :: ts
     logical                 :: startFlag
     character*100           :: odir
     integer                 :: nsites
     integer                 :: ntimes
     real,   allocatable     :: time(:)
     real,   allocatable     :: WLobs(:,:)
     real,   allocatable     :: sites_data(:,:)
  end type hydrowebwlobsdec

  type(hydrowebwlobsdec), allocatable:: hydrowebWLobs(:)

contains
  
!BOP
! 
! !ROUTINE: hydrowebWL_obsInit
! \label{hydrowebWL_obsInit}
! 
! !INTERFACE: 
  subroutine hydrowebWL_obsinit()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_timeMgrMod
    use LDT_logMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading the LPRM vegetation optical depth data. 
! 
!EOP
    integer                 :: npts
    integer                 :: ftn
    integer                 :: siteid
    character*100           :: wldistancemap
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(hydrowebWLobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'Hydroweb water level data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, hydrowebWLobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'Hydroweb water level data directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Hydroweb water level distance map:",&
         rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,wldistancemap,&
            rc=status)
       call LDT_verify(status, &
            'Hydroweb water level data distance map: not defined')
    enddo

    do n=1,LDT_rc%nnest
       call ESMF_TimeIntervalSet(hydrowebwlobs(n)%ts,d=1,rc=status)
       call LDT_verify(status, 'Error in ESMF_TimeIntervalSet in hydrowebwlobs')

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%wl_obs, &
            "m",1,1)
       LDT_DAobsData(n)%wl_obs%selectStats = 1
    
       allocate(hydrowebWLobs(n)%sites_data(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       hydrowebWLobs(n)%sites_data = -9999.0


#if(defined USE_NETCDF3 || defined USE_NETCDF4)
       call LDT_verify(nf90_open(path=trim(wldistancemap),&
            mode=NF90_NOWRITE,ncid=ftn),&
            'Error opening file '//trim(wldistancemap))
       call LDT_verify(nf90_inq_varid(ftn,'sites',siteid),&
            'Error with nf90_inq_varid: sites')
       call LDT_verify(nf90_get_var(ftn,siteid,&
            hydrowebWLobs(n)%sites_data, &
            start=(/LDT_ews_halo_ind(n,LDT_localPet+1),&
            LDT_nss_halo_ind(n,LDT_localPet+1)/),&
            count = (/LDT_rc%lnc(n),LDT_rc%lnr(n)/)),&
            'Error with nf90_get_var: sites')
       call LDT_verify(nf90_close(ftn))
#endif

       hydrowebWLobs(n)%startFlag = .true.

       call ESMF_TimeIntervalSet(hydrowebwlobs(n)%ts,d=1,rc=status)
       call LDT_verify(status, 'Error in ESMF_TimeIntervalSet in hydrowebwlobs')
 
    enddo
  end subroutine hydrowebWL_obsinit
     
end module hydrowebWL_obsMod
