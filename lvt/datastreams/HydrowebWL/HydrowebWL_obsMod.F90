!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: HydrowebWL_obsMod
!  \label(HydrowebWL_obsMod)
!
! !INTERFACE:
module HydrowebWL_obsMod
! 
! !USES: 
  use ESMF

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the HydrowebWL
!  in-situ data. 
!  
!  web link: http://www.ozflux.org.au
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Apr 2020;   Sujay Kumar  Initial Specification
! 
!EOP
! 
! 
!
  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: HydrowebWL_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: HydrowebWLobs
!EOP
  type, public :: hydrowebwlobsdec
     integer                     :: nsites
     integer                     :: da
     real, allocatable           :: sites(:,:)
     real, allocatable           :: wl(:,:)
     type(ESMF_Time)             :: startTime
     type(ESMF_TimeInterval)     :: ts

  end type hydrowebwlobsdec

  type(hydrowebwlobsdec), allocatable:: hydrowebwlobs(:)

contains
  
!BOP
! 
! !ROUTINE: HydrowebWL_obsInit
! \label{HydrowebWL_obsInit}
!
! !INTERFACE: 
  subroutine HydrowebWL_obsinit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod
    use map_utils
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading HydrowebWL data. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!! !ARGUMENTS: 
    integer,   intent(IN) :: i 
!!EOP
    integer                 :: status
    integer                 :: ftn
    integer                 :: ncId, nrId, tId, wlid,siteid
    integer                 :: nc,nr,nt
    character*100           :: stnlist_file
    character*100           :: odir
    integer                 :: c,r,latid, lonid
    real,    allocatable    :: lat(:,:),lon(:,:)

    if(.not.allocated(hydrowebwlobs)) then 
       allocate(hydrowebwlobs(LVT_rc%nDataStreams))
    endif
!

    call ESMF_ConfigGetAttribute(LVT_Config, odir, &
         label='Hydroweb observation directory:', rc=status)
    call LVT_verify(status, 'Hydroweb observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, stnlist_file, &
         label='Hydroweb station list file:', rc=status)
    call LVT_verify(status, 'Hydroweb station list file: not defined')

    write(LVT_logunit,*) '[INFO] Processing Hydroweb data locations '

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    status = nf90_open(path=trim(stnlist_file),mode=NF90_NOWRITE,ncid=ftn)
    call LVT_verify(status, 'Error opening file'//trim(stnlist_file))

    status = nf90_inq_dimid(ftn, 'east_west',ncId)
    call LVT_verify(status, 'Error nf90_inq_dimid: east_west')
    
    status = nf90_inquire_dimension(ftn, ncId, len=nc)
    call LVT_verify(status, 'Error nf90_inquire_dimension: east_west')

    status = nf90_inq_dimid(ftn, 'north_south',nrId)
    call LVT_verify(status, 'Error nf90_inq_dimid: north_south')
    
    status = nf90_inquire_dimension(ftn, nrId, len=nr)
    call LVT_verify(status, 'Error nf90_inquire_dimension: north_south')

    allocate(hydrowebwlobs(i)%sites(nc,nr))
    hydrowebwlobs(i)%sites = -1

    allocate(lat(nc,nr))
    allocate(lon(nc,nr))
    !variable ids
    status = nf90_inq_varid(ftn,'lat',latid)
    call LVT_verify(status,'Error in nf90_inq_varid: lat')

    status = nf90_inq_varid(ftn,'lon',lonid)
    call LVT_verify(status,'Error in nf90_inq_varid: lon')

    status = nf90_get_var(ftn,latid, lat)
    call LVT_verify(status, 'Error nf90_get_var: lat')

    status = nf90_get_var(ftn,lonid, lon)
    call LVT_verify(status, 'Error nf90_get_var: lon')

    status = nf90_inq_varid(ftn,'sites',siteid)
    call LVT_verify(status,'Error in nf90_inq_varid: sites')
    
    status = nf90_get_var(ftn,siteid, hydrowebwlobs(i)%sites)
    call LVT_verify(status, 'Error nf90_get_var: sites')

    status = nf90_close(ftn)
    call LVT_verify(status,'Error in nf90_close')

#endif

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    status = nf90_open(path=trim(odir),mode=NF90_NOWRITE,ncid=ftn)
    call LVT_verify(status, 'Error opening file'//trim(odir))

    status = nf90_inq_dimid(ftn, 'time',tId)
    call LVT_verify(status, 'Error nf90_inq_dimid: time')
    
    status = nf90_inquire_dimension(ftn, tId, len=nt)
    call LVT_verify(status, 'Error nf90_inquire_dimension: time')

    status = nf90_inq_dimid(ftn, 'sites',siteId)
    call LVT_verify(status, 'Error nf90_inq_dimid: sites')
    
    status = nf90_inquire_dimension(ftn, siteId, len=hydrowebwlobs(i)%nsites)
    call LVT_verify(status, 'Error nf90_inquire_dimension: sites')

    allocate(hydrowebwlobs(i)%wl(hydrowebwlobs(i)%nsites,nt))

    !variable ids
    status = nf90_inq_varid(ftn,'water_elevation',wlid)
    call LVT_verify(status,'Error in nf90_inq_varid: water_elevation')

    status = nf90_get_var(ftn,wlid, hydrowebwlobs(i)%wl)
    call LVT_verify(status, 'Error nf90_get_var: water_elevation')

    status = nf90_close(ftn)
    call LVT_verify(status,'Error in nf90_close')
#endif

     call ESMF_TimeSet(hydrowebwlobs(i)%startTime, yy=2016,&
          mm = 7,&
          dd = 1,&
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: HydrowebWLobs')

     call ESMF_TimeIntervalSet(hydrowebwlobs(i)%ts,d=1,rc=status)
     call LVT_verify(status, 'error in timeintervalset: HydrowebWLobs')

     hydrowebwlobs(i)%da = -1

    call LVT_update_timestep(LVT_rc, 86400)
  end subroutine HydrowebWL_obsinit

end module HydrowebWL_obsMod
