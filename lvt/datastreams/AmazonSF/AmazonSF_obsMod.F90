!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: AmazonSF_obsMod
!  \label(AmazonSF_obsMod)
!
! !INTERFACE:
module AmazonSF_obsMod
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
!  This module handles the observation plugin for the AmazonSF
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
  PUBLIC :: AmazonSF_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: AmazonSFobs
!EOP
  type, public :: amazonsfobsdec
     integer                     :: nsites
     integer                     :: da
     real, allocatable           :: sites(:,:)
     real, allocatable           :: sf(:,:)
     type(ESMF_Time)             :: startTime
     type(ESMF_TimeInterval)     :: ts

  end type amazonsfobsdec

  type(amazonsfobsdec), allocatable:: amazonsfobs(:)

contains
  
!BOP
! 
! !ROUTINE: AmazonSF_obsInit
! \label{AmazonSF_obsInit}
!
! !INTERFACE: 
  subroutine AmazonSF_obsinit(i)
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
!  for reading AmazonSF data. 
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
    integer                 :: ncId, nrId, tId, sfid,siteid
    integer                 :: nc,nr,nt
    character*100           :: stnlist_file
    character*100           :: odir
    integer                 :: c,r,latid, lonid
    real,    allocatable    :: lat(:,:),lon(:,:)

    if(.not.allocated(amazonsfobs)) then 
       allocate(amazonsfobs(LVT_rc%nDataStreams))
    endif
!

    call ESMF_ConfigGetAttribute(LVT_Config, odir, &
         label='AmazonSF observation directory:', rc=status)
    call LVT_verify(status, 'AmazonSF observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, stnlist_file, &
         label='AmazonSF station list file:', rc=status)
    call LVT_verify(status, 'AmazonSF station list file: not defined')

    write(LVT_logunit,*) '[INFO] Processing AmazonSF data locations '

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

    allocate(amazonsfobs(i)%sites(nc,nr))
    amazonsfobs(i)%sites = -1

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
    
    status = nf90_get_var(ftn,siteid, amazonsfobs(i)%sites)
    call LVT_verify(status, 'Error nf90_get_var: sites')

    status = nf90_close(ftn)
    call LVT_verify(status,'Error in nf90_close')

!    do r=1,nr
!       do c=1,nc
!          if(amazonsfobs(i)%sites(c,r).gt.0) then
!             print*, lat(c,r), lon(c,r), lat(c,r), lon(c,r)
!          endif
!       enddo
!    enddo
!    stop
#if 0 
    c=102
    r=81
    print*, c,r,lat(c,r),lon(c,r)
    c=72
    r=72
    print*, c,r,lat(c,r),lon(c,r)
    c=82
    r=75
    print*, c,r,lat(c,r),lon(c,r)
    c=91
    r=76
    print*, c,r,lat(c,r),lon(c,r)
    c=47
    r=83
    print*, c,r,lat(c,r),lon(c,r)
    stop
#endif
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
    
    status = nf90_inquire_dimension(ftn, siteId, len=amazonsfobs(i)%nsites)
    call LVT_verify(status, 'Error nf90_inquire_dimension: sites')

    allocate(amazonsfobs(i)%sf(amazonsfobs(i)%nsites,nt))

    !variable ids
    status = nf90_inq_varid(ftn,'streamflow',sfid)
    call LVT_verify(status,'Error in nf90_inq_varid: streamflow')

    status = nf90_get_var(ftn,sfid, amazonsfobs(i)%sf)
    call LVT_verify(status, 'Error nf90_get_var: streamflow')

    status = nf90_close(ftn)
    call LVT_verify(status,'Error in nf90_close')
#endif

     call ESMF_TimeSet(amazonsfobs(i)%startTime, yy=2016,&
          mm = 1,&
          dd = 12,&
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: AmazonSFobs')

     call ESMF_TimeIntervalSet(amazonsfobs(i)%ts,d=1,rc=status)
     call LVT_verify(status, 'error in timeintervalset: AmazonSFobs')

     amazonsfobs(i)%da = -1

    call LVT_update_timestep(LVT_rc, 86400)
  end subroutine AmazonSF_obsinit

end module AmazonSF_obsMod
