!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: MCD15A2Hlai_Mod
! 
! !DESCRIPTION: 
!   This module contains interfaces and subroutines to
!   handle the MODIS MCD15A2H LAI data.
! 
! !REVISION HISTORY: 
!  16 Jun 2020    Wanshu Nie; initial specification
! 
module MCD15A2Hlai_Mod
! !USES: 
  use ESMF
  use map_utils
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: MCD15A2Hlai_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: MCD15A2Hlai_struc
!EOP
  type, public:: MCD15A2Hlai_dec

     character*100          :: version
     logical                :: startMode
     integer                :: nc
     integer                :: nr
     integer                :: mi
     real,     allocatable  :: laiobs1(:)
     real,     allocatable  :: laiobs2(:)
     real                   :: gridDesci(50)    
     real*8                 :: time1, time2
     integer                :: fnd
     integer                :: tsmooth
     integer                :: climofill
     integer                :: qcflag
     real,    allocatable :: rlat(:)
     real,    allocatable :: rlon(:)
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type MCD15A2Hlai_dec
  
  type(MCD15A2Hlai_dec),allocatable :: MCD15A2Hlai_struc(:)
  
contains

!BOP
! 
! !ROUTINE: MCD15A2Hlai_setup
! \label{MCD15A2Hlai_setup}
! 
! !INTERFACE: 
  subroutine MCD15A2Hlai_setup(k, OBS_State, OBS_Pert_State)
! !USES: 
    use LIS_coreMod
    use LIS_timeMgrMod
    use LIS_historyMod
    use LIS_dataAssimMod
    use LIS_perturbMod
    use LIS_DAobservationsMod
    use LIS_logmod

    implicit none 

! !ARGUMENTS: 
    integer                ::  k
    type(ESMF_State)       ::  OBS_State(LIS_rc%nnest)
    type(ESMF_State)       ::  OBS_Pert_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This routine completes the runtime initializations and 
!   creation of data structures required for handling MCD15A2H LAI data.
!  
!   The arguments are: 
!   \begin{description}
!    \item[k] number of observation state 
!    \item[OBS\_State]   observation state 
!    \item[OBS\_Pert\_State] observation perturbations state
!   \end{description}
!EOP
    integer                ::  n,i,t,kk,jj
    integer                ::  ftn
    integer                ::  status
    type(ESMF_Field)       ::  obsField(LIS_rc%nnest)
    type(ESMF_Field)       ::  pertField(LIS_rc%nnest)
    type(ESMF_ArraySpec)   ::  intarrspec, realarrspec
    type(ESMF_ArraySpec)   ::  pertArrSpec
    character(len=LIS_CONST_PATH_LEN) ::  laiobsdir
    character*100          ::  temp
    real,  allocatable         ::  ssdev(:)
    character*1            ::  vid(2)
    type(pert_dec_type)    ::  obs_pert
    real, pointer          ::  obs_temp(:,:)
    character*40, allocatable  ::  vname(:)
    real        , allocatable  ::  varmin(:)
    real        , allocatable  ::  varmax(:)
    integer                :: c,r
    real*8                 ::  cornerlat1, cornerlon1
    real*8                 ::  cornerlat2, cornerlon2

    allocate(MCD15A2Hlai_struc(LIS_rc%nnest))

    call ESMF_ArraySpecSet(intarrspec,rank=1,typekind=ESMF_TYPEKIND_I4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(realarrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ArraySpecSet(pertArrSpec,rank=2,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI data directory:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,laiobsdir,&
            rc=status)
       call LIS_verify(status, 'MCD15A2H LAI data directory: is missing')

       call ESMF_AttributeSet(OBS_State(n),"Data Directory",&
            laiobsdir, rc=status)
       call LIS_verify(status)
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI data version:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MCD15A2Hlai_struc(n)%version,&
            rc=status)
       call LIS_verify(status, 'MCD15A2H LAI data version: is missing')
    enddo


    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI apply temporal smoother between 8-day intervals:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MCD15A2Hlai_struc(n)%tsmooth,&
            rc=status)
       call LIS_verify(status, 'MCD15A2H LAI apply temporal smoother between 8-day intervals: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI apply climatological fill values:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MCD15A2Hlai_struc(n)%climofill,&
            rc=status)
       call LIS_verify(status, 'MCD15A2H LAI apply climatological fill values: is missing')
    enddo

    call ESMF_ConfigFindLabel(LIS_config,"MCD15A2H LAI apply QC flags:",&
         rc=status)
    do n=1,LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MCD15A2Hlai_struc(n)%qcflag,&
            rc=status)
       call LIS_verify(status, 'MCD15A2H LAI apply QC flags: is missing')
    enddo

   do n=1,LIS_rc%nnest
       call ESMF_AttributeSet(OBS_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Update Time",&
            -99.0, rc=status)
       call LIS_verify(status)

       call ESMF_AttributeSet(OBS_State(n),"Data Assimilate Status",&
            .false., rc=status)
       call LIS_verify(status)
       
       call ESMF_AttributeSet(OBS_State(n),"Number Of Observations",&
            LIS_rc%obs_ngrid(k),rc=status)
       call LIS_verify(status)
       
    enddo

    write(LIS_logunit,*)&
         '[INFO] read MCD15A2H LAI data specifications'       

    do n=1,LIS_rc%nnest

       allocate(MCD15A2Hlai_struc(n)%laiobs1(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       allocate(MCD15A2Hlai_struc(n)%laiobs2(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
       
       write(unit=temp,fmt='(i2.2)') 1
       read(unit=temp,fmt='(2a1)') vid

       obsField(n) = ESMF_FieldCreate(arrayspec=realarrspec,&
            grid=LIS_obsVecGrid(n,k),&
            name="Observation"//vid(1)//vid(2), rc=status)
       call LIS_verify(status)

!Perturbations State
       write(LIS_logunit,*) '[INFO] Opening attributes for observations ',&
            trim(LIS_rc%obsattribfile(k))
       ftn = LIS_getNextUnitNumber()
       open(ftn,file=trim(LIS_rc%obsattribfile(k)),status='old')
       read(ftn,*)
       read(ftn,*) LIS_rc%nobtypes(k)
       read(ftn,*)
    
       allocate(vname(LIS_rc%nobtypes(k)))
       allocate(varmax(LIS_rc%nobtypes(k)))
       allocate(varmin(LIS_rc%nobtypes(k)))
       
       do i=1,LIS_rc%nobtypes(k)
          read(ftn,fmt='(a40)') vname(i)
          read(ftn,*) varmin(i),varmax(i)
          write(LIS_logunit,*) '[INFO] ',vname(i),varmin(i),varmax(i)
       enddo
       call LIS_releaseUnitNumber(ftn)   

       allocate(ssdev(LIS_rc%obs_ngrid(k)))
      
       if(trim(LIS_rc%perturb_obs(k)).ne."none") then 
          allocate(obs_pert%vname(1))
          allocate(obs_pert%perttype(1))
          allocate(obs_pert%ssdev(1))
          allocate(obs_pert%stdmax(1))
          allocate(obs_pert%zeromean(1))
          allocate(obs_pert%tcorr(1))
          allocate(obs_pert%xcorr(1))
          allocate(obs_pert%ycorr(1))
          allocate(obs_pert%ccorr(1,1))
          
          call LIS_readPertAttributes(1,LIS_rc%obspertAttribfile(k),&
               obs_pert)
          
          ! Set obs err to be uniform (will be rescaled later for each grid point). 
          ssdev = obs_pert%ssdev(1)
          
          pertField(n) = ESMF_FieldCreate(arrayspec=pertArrSpec,&
               grid=LIS_obsEnsOnGrid(n,k),name="Observation"//vid(1)//vid(2),&
               rc=status)
          call LIS_verify(status)
          
          ! initializing the perturbations to be zero 
          call ESMF_FieldGet(pertField(n),localDE=0,farrayPtr=obs_temp,rc=status)
          call LIS_verify(status)
          obs_temp(:,:) = 0 
          
          call ESMF_AttributeSet(pertField(n),"Perturbation Type",&
               obs_pert%perttype(1), rc=status)
          call LIS_verify(status)
          
          if(LIS_rc%obs_ngrid(k).gt.0) then 
             call ESMF_AttributeSet(pertField(n),"Standard Deviation",&
                  ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
             call LIS_verify(status)
          endif
          
          call ESMF_AttributeSet(pertField(n),"Std Normal Max",&
               obs_pert%stdmax(1), rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Ensure Zero Mean",&
               obs_pert%zeromean(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"Temporal Correlation Scale",&
               obs_pert%tcorr(1),rc=status)
          call LIS_verify(status)
          
          call ESMF_AttributeSet(pertField(n),"X Correlation Scale",&
               obs_pert%xcorr(1),rc=status)
          
          call ESMF_AttributeSet(pertField(n),"Y Correlation Scale",&
               obs_pert%ycorr(1),rc=status)

          call ESMF_AttributeSet(pertField(n),"Cross Correlation Strength",&
               obs_pert%ccorr(1,:),itemCount=1,rc=status)

          call ESMF_StateAdd(OBS_Pert_State(n),(/pertField(n)/),rc=status)
          call LIS_verify(status)         

       endif
       
       deallocate(vname)
       deallocate(varmax)
       deallocate(varmin)
       deallocate(ssdev)   
       
    enddo
    
    do n=1,LIS_rc%nnest

       if(LIS_rc%lis_obs_map_proj(k).eq."latlon") then
          cornerlat1 = max(-59.9978927, nint((LIS_rc%obs_gridDesc(k,4)+59.9978927)/0.00416667)*0.00416667-59.9978927-50*0.00416667)
          cornerlon1 = max(-179.9979167, nint((LIS_rc%obs_gridDesc(k,5)+179.9979167)/0.00416667)*0.00416667-179.9979167-50*0.00416667)
          cornerlat2 = min(89.9979167, nint((LIS_rc%obs_gridDesc(k,7)+59.9978927)/0.00416667)*0.00416667-59.9978927+50*0.00416667)
          cornerlon2 = min(179.9979167, nint((LIS_rc%obs_gridDesc(k,8)+179.9979167)/0.00416667)*0.00416667-179.9979167+50*0.00416667)
       elseif(LIS_rc%lis_obs_map_proj(k).eq."lambert") then
          cornerlat1 = max(-59.9978927, nint((LIS_rc%minLat(n)+59.9978927)/0.00416667)*0.00416667-59.9978927-50*0.00416667)
          cornerlon1 = max(-179.9979167, nint((LIS_rc%minLon(n)+179.9979167)/0.00416667)*0.00416667-179.9979167-50*0.00416667)
          cornerlat2 = min(89.9979167, nint((LIS_rc%maxLat(n)+59.9978927)/0.00416667)*0.00416667-59.9978927+50*0.00416667)
          cornerlon2 = min(179.9979167, nint((LIS_rc%maxLon(n)+179.9979167)/0.00416667)*0.00416667-179.9979167+50*0.00416667)
       endif

       MCD15A2Hlai_struc(n)%nc = nint((cornerlon2-cornerlon1)/0.00416667)+1
       MCD15A2Hlai_struc(n)%nr = nint((cornerlat2-cornerlat1)/0.00416667)+1

       MCD15A2Hlai_struc(n)%gridDesci(1) = 0 
       MCD15A2Hlai_struc(n)%gridDesci(2) = MCD15A2Hlai_struc(n)%nc
       MCD15A2Hlai_struc(n)%gridDesci(3) = MCD15A2Hlai_struc(n)%nr 
       MCD15A2Hlai_struc(n)%gridDesci(4) = cornerlat1
       MCD15A2Hlai_struc(n)%gridDesci(5) = cornerlon1
       MCD15A2Hlai_struc(n)%gridDesci(6) = 128
       MCD15A2Hlai_struc(n)%gridDesci(7) = cornerlat2
       MCD15A2Hlai_struc(n)%gridDesci(8) = cornerlon2
       MCD15A2Hlai_struc(n)%gridDesci(9) = 0.00416667
       MCD15A2Hlai_struc(n)%gridDesci(10) = 0.00416667
       MCD15A2Hlai_struc(n)%gridDesci(20) = 64

       MCD15A2Hlai_struc(n)%mi = MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr

!-----------------------------------------------------------------------------
!   Use interpolation if LIS is running finer than 500 m. 
!-----------------------------------------------------------------------------
       if(LIS_rc%obs_gridDesc(k,10).le.0.00416667) then 

          allocate(MCD15A2Hlai_struc(n)%rlat(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%rlon(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%n11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%n12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%n21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%n22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%w11(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%w12(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%w21(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          allocate(MCD15A2Hlai_struc(n)%w22(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k)))
          
          call bilinear_interp_input_withgrid(MCD15A2Hlai_struc(n)%gridDesci(:), &
               LIS_rc%obs_gridDesc(k,:),&
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
               MCD15A2Hlai_struc(n)%rlat, MCD15A2Hlai_struc(n)%rlon,&
               MCD15A2Hlai_struc(n)%n11, MCD15A2Hlai_struc(n)%n12, &
               MCD15A2Hlai_struc(n)%n21, MCD15A2Hlai_struc(n)%n22, &
               MCD15A2Hlai_struc(n)%w11, MCD15A2Hlai_struc(n)%w12, &
               MCD15A2Hlai_struc(n)%w21, MCD15A2Hlai_struc(n)%w22)
       else
          
          allocate(MCD15A2Hlai_struc(n)%n11(&
               MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr))

          call upscaleByAveraging_input(MCD15A2Hlai_struc(n)%gridDesci(:),&
               LIS_rc%obs_gridDesc(k,:),&
               MCD15A2Hlai_struc(n)%nc*MCD15A2Hlai_struc(n)%nr, &
               LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), MCD15A2Hlai_struc(n)%n11)
       endif
       
       call LIS_registerAlarm("MCD15A2H LAI read alarm",&
            86400.0, 86400.0)

       MCD15A2Hlai_struc(n)%startMode = .true. 

       call ESMF_StateAdd(OBS_State(n),(/obsField(n)/),rc=status)
       call LIS_verify(status)
     
    enddo
  end subroutine MCD15A2Hlai_setup
end module MCD15A2Hlai_Mod
