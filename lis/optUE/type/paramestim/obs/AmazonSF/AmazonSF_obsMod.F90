!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !MODULE: AmazonSF_obsMod
! 
! !DESCRIPTION: 
!  This module handles the observation plugin for the
!  streamflow data over Amazon for use 
!  within the LIS OPT/UE subsystem. 
!
!   
! !REVISION HISTORY: 
!  2 Oct 2020 Sujay Kumar;   Initial Specification 
! 
module AmazonSF_obsMod
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: AmazonSF_obs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: AmazonSF_obs_struc

  type, public ::  AmazonSF_obs_data_dec

     integer                    :: nsites
     integer                    :: da
     type(ESMF_Time)            :: starttime
     type(ESMF_TimeInterval)    :: timestep
     
     real, allocatable          :: sites(:,:)
     real, allocatable          :: sf(:,:)

  end type AmazonSF_obs_data_dec

  type(AmazonSF_obs_data_dec), allocatable :: AmazonSF_obs_struc(:)

contains
!BOP
! 
! !ROUTINE: AmazonSF_obs_setup
! \label{AmazonSF_obs_setup}
! 
! !INTERFACE: 
  subroutine AmazonSF_obs_setup(Obs_State)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!  This routines completes the setup of the Amazon streamflow plugin
!  for OPTUE. This includes the definition of the data grid and the 
!  setup of the interpolation weights. 
!   
!   The arguments are: 
!   \begin{description}
!    \item[Obs\_State]   observation state object 
!   \end{description}
!EOP
    integer                   :: n 
    integer                   :: ncId, nrId,sfid
    integer                   :: latId,lonId,siteid,tid
    integer                   :: nc,nr,nt
    type(ESMF_ArraySpec)      :: realarrspec
    type(ESMF_Field)          :: obsField
    character*100             :: obsdir,stnlist_file
    real,  allocatable        :: lat(:,:),lon(:,:),sites(:,:)
    integer                   :: k
    integer                   :: ftn
    integer                   :: status


    allocate(AmazonSF_obs_struc(LIS_rc%nnest))

    write(LIS_logunit,*) '[INFO] Setting up Amazon streamflow data reader....'

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigFindLabel(LIS_config,"Amazon streamflow data directory:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,obsdir,&
         rc=status)
    call LIS_verify(status, 'Amazon streamflow data directory: not defined')
    
    call ESMF_ConfigFindLabel(LIS_config,"Amazon streamflow station list file:",&
         rc=status)
    call ESMF_ConfigGetAttribute(LIS_config,stnlist_file,&
         rc=status)
    call LIS_verify(status, 'Amazon streamflow station list file: not defined')
    
    do n=1,LIS_rc%nnest

       call ESMF_AttributeSet(Obs_State(n),"Data Directory",&
            obsdir, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

    do n=1,LIS_rc%nnest

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

       status = nf90_open(path=trim(stnlist_file),mode=NF90_NOWRITE,ncid=ftn)
       call LIS_verify(status, 'Error opening file'//trim(stnlist_file))
       
       status = nf90_inq_dimid(ftn, 'east_west',ncId)
       call LIS_verify(status, 'Error nf90_inq_dimid: east_west')
       
       status = nf90_inquire_dimension(ftn, ncId, len=nc)
       call LIS_verify(status, 'Error nf90_inquire_dimension: east_west')
       
       status = nf90_inq_dimid(ftn, 'north_south',nrId)
       call LIS_verify(status, 'Error nf90_inq_dimid: north_south')
       
       status = nf90_inquire_dimension(ftn, nrId, len=nr)
       call LIS_verify(status, 'Error nf90_inquire_dimension: north_south')
       
       allocate(sites(nc,nr))
       sites = -1
       
       allocate(lat(nc,nr))
       allocate(lon(nc,nr))
       !variable ids
       status = nf90_inq_varid(ftn,'lat',latid)
       call LIS_verify(status,'Error in nf90_inq_varid: lat')
       
       status = nf90_inq_varid(ftn,'lon',lonid)
       call LIS_verify(status,'Error in nf90_inq_varid: lon')
       
       status = nf90_get_var(ftn,latid, lat)
       call LIS_verify(status, 'Error nf90_get_var: lat')
       
       status = nf90_get_var(ftn,lonid, lon)
       call LIS_verify(status, 'Error nf90_get_var: lon')
       
       status = nf90_inq_varid(ftn,'sites',siteid)
       call LIS_verify(status,'Error in nf90_inq_varid: sites')
       
       status = nf90_get_var(ftn,siteid, sites)
       call LIS_verify(status, 'Error nf90_get_var: sites')
       
       status = nf90_close(ftn)
       call LIS_verify(status,'Error in nf90_close')

       allocate(AmazonSF_obs_struc(n)%sites(LIS_rc%lnc(n),LIS_rc%lnr(n)))
       
       AmazonSF_obs_struc(n)%sites(:,:) = sites(&
          LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))

       deallocate(sites)


       status = nf90_open(path=trim(obsdir),mode=NF90_NOWRITE,ncid=ftn)
       call LIS_verify(status, 'Error opening file'//trim(obsdir))
       
       status = nf90_inq_dimid(ftn, 'time',tId)
       call LIS_verify(status, 'Error nf90_inq_dimid: time')
       
       status = nf90_inquire_dimension(ftn, tId, len=nt)
       call LIS_verify(status, 'Error nf90_inquire_dimension: time')
       
       status = nf90_inq_dimid(ftn, 'sites',siteId)
       call LIS_verify(status, 'Error nf90_inq_dimid: sites')
       
       status = nf90_inquire_dimension(ftn, siteId, len=AmazonSF_obs_struc(n)%nsites)
       call LIS_verify(status, 'Error nf90_inquire_dimension: sites')
       
       allocate(AmazonSF_obs_struc(n)%sf(AmazonSF_obs_struc(n)%nsites,nt))
       
       !variable ids
       status = nf90_inq_varid(ftn,'streamflow',sfid)
       call LIS_verify(status,'Error in nf90_inq_varid: streamflow')
       
       status = nf90_get_var(ftn,sfid, AmazonSF_obs_struc(n)%sf)
       call LIS_verify(status, 'Error nf90_get_var: streamflow')
       
       status = nf90_close(ftn)
       call LIS_verify(status,'Error in nf90_close')
       
#endif

       call ESMF_TimeSet(AmazonSF_obs_struc(n)%startTime, yy=2016,&
          mm = 1,&
          dd = 12,&
          h = 0, &
          m = 0, &
          calendar = LIS_calendar, &
          rc=status)
       call LIS_verify(status,'error in timeset: AmazonSFobs')

       obsField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecRoutingGrid(n), &
            name="Amazon_SF", rc=status)
       call LIS_verify(status)
       
       call ESMF_StateAdd(Obs_State(n),(/obsField/),rc=status)
       call LIS_verify(status)

       call ESMF_TimeIntervalSet(AmazonSF_obs_struc(n)%timestep, &
            s=86400,rc=status)
       call LIS_verify(status, 'Error in ESMF_TimeIntervalSet (Amazon SFobs)')
       
       AmazonSF_obs_struc(n)%da = -1
    enddo

    write(LIS_logunit,*) '[INFO] created the States to hold the Amazon streamflow data'
    
  end subroutine AmazonSF_obs_setup
  
end module AmazonSF_obsMod
