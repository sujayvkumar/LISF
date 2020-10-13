!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module HYMAP2_peMod
!BOP
!
! !MODULE: HYMAP2_peMod
!
! !DESCRIPTION:
!  This module contains the definitions of the HYMAP2 model parameters
!  used in parameter estimation. The data structure is used to expose
!  the ROUTING parameters to be used in opt/ue. 
!
! !REVISION HISTORY:
!  2 Oct 2020; Sujay Kumar, Initial Code
! !USES:        
  use ESMF
  use LIS_numerRecipesMod, only : LIS_rand_func

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: HYMAP2_setup_pedecvars
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: HYMAP2_pe_struc

!EOP
  type, public ::  HYMAP2_pe_dec 
     integer               :: nparams
     character*40, allocatable :: param_name(:)
     integer     , allocatable :: param_select(:)
     real        , allocatable :: param_min(:)
     real        , allocatable :: param_max(:)
  end type HYMAP2_pe_dec

  type(HYMAP2_pe_dec), allocatable :: HYMAP2_pe_struc(:)

  SAVE
contains

!BOP
! !ROUTINE: HYMAP2_setup_pedecvars
!  \label{HYMAP2_setup_pedecvars}
!
! !REVISION HISTORY:
!  2 Oct 2020; Sujay Kumar, Initial Code
!
! !INTERFACE:
  subroutine HYMAP2_setup_pedecvars(DEC_State, Feas_State)
! !USES:
    use ESMF
    use LIS_coreMod
    use LIS_logMod
    use HYMAP2_routingMod,     only : HYMAP2_routing_struc

    implicit none
! !ARGUMENTS: 
    character*100               :: decSpaceAttribsFile
    type(ESMF_State)            :: DEC_State
    type(ESMF_State)            :: Feas_State

!
! !DESCRIPTION:
!  
!  This routine determines the list of parameters to be used in parameter
!  estimation, initializes them, and updates the LIS decision space.  
! 
!EOP

    integer                     :: n 
    type(ESMF_ArraySpec)        :: arrspec1
    type(ESMF_Field)            :: varField
    type(ESMF_Field)            :: feasField
    real, pointer               :: vardata(:)
    integer, pointer            :: mod_flag(:)
    integer                     :: i,t,k,m
    integer                     :: status
    integer, parameter          :: seed_base=-1000
    integer                     :: seed
    real                        :: rand
    integer                     :: NT    
    character*100               :: vname
    integer                     :: count
    integer                     :: gid

    call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
    call LIS_verify(status)
    call ESMF_FieldGet(feasField,localDE=0,farrayPtr=mod_flag,rc=status)
    call LIS_verify(status)
    
    call ESMF_ConfigGetAttribute(LIS_config,decSpaceAttribsFile,&
         label="Routing decision space attributes file:",rc=status)
    call LIS_verify(status, "Routing decision space attributes file: not defined")

    allocate(HYMAP2_pe_struc(LIS_rc%nnest))
    n = 1
    HYMAP2_pe_struc(n)%nparams = 3

    allocate(HYMAP2_pe_struc(n)%param_name(HYMAP2_pe_struc(n)%nparams))
    allocate(HYMAP2_pe_struc(n)%param_select(HYMAP2_pe_struc(n)%nparams))
    allocate(HYMAP2_pe_struc(n)%param_min(HYMAP2_pe_struc(n)%nparams))
    allocate(HYMAP2_pe_struc(n)%param_max(HYMAP2_pe_struc(n)%nparams))

    ! read the attributes file. 
    call LIS_readPEDecSpaceAttributes(decSpaceAttribsFile, &
         HYMAP2_pe_struc(n)%nparams, &
         HYMAP2_pe_struc(n)%param_name, &
         HYMAP2_pe_struc(n)%param_select, &
         HYMAP2_pe_struc(n)%param_min, &
         HYMAP2_pe_struc(n)%param_max)

    call ESMF_ArraySpecSet(arrspec1,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)
    
    do i=1,HYMAP2_pe_struc(n)%nparams
       if(HYMAP2_pe_struc(n)%param_select(i).eq.1) then 
          vname=trim(HYMAP2_pe_struc(n)%param_name(i))
          varField = ESMF_FieldCreate(arrayspec=arrspec1, &
               grid=LIS_vecRoutingTile(n),&
               name=vname,&
               rc=status)
          call LIS_verify(status, &
               'problem with fieldcreate in HYMAP2_setup_pedecvars')
          call ESMF_AttributeSet(varField,'MinRange',&
               HYMAP2_pe_struc(n)%param_min(i),rc=status)
          call LIS_verify(status, &
               'setting minrange to decspace obj in HYMAP2_setup_devars')
          call ESMF_AttributeSet(varField,'MaxRange',&
               HYMAP2_pe_struc(n)%param_max(i),rc=status)
          call LIS_verify(status, &
               'setting maxrange to decspace obj in HYMAP2_setup_devars')

          call ESMF_StateAdd(DEC_State,(/varField/),rc=status)
          call LIS_verify(status,&
               'stateadd in HYMAP2_setup_pedecvars')

          call ESMF_StateGet(DEC_State,vname,varField,rc=status)
          call LIS_verify(status)
          
          call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
               rc=status)
          call LIS_verify(status)
          

          NT=HYMAP2_routing_struc(n)%nseqall
          if(vname.eq."RIVWTH")  then  
             do t=1,NT
                do m=1,LIS_rc%nensem(n)
                   k=(t-1)*LIS_rc%nensem(n)+m
                   vardata(k) = HYMAP2_routing_struc(n)%rivwth(t,m)
                enddo
             enddo
          endif

          if(vname.eq."RIVHGT")  then 
             do t=1,NT
                do m=1,LIS_rc%nensem(n)
                   k=(t-1)*LIS_rc%nensem(n)+m
                   vardata(k) = HYMAP2_routing_struc(n)%rivhgt(t,m)
                enddo
             enddo
          endif
          if(vname.eq."RIVROUGHNESS") then 
             do t=1,NT
                do m=1,LIS_rc%nensem(n)
                   k=(t-1)*LIS_rc%nensem(n)+m
                   vardata(k) = HYMAP2_routing_struc(n)%rivman(t,m)
                enddo
             enddo
          endif

          !Test whether any defaults are out of bounds
          count=0
          do t=1,HYMAP2_routing_struc(n)%nseqall/LIS_rc%nensem(n)
             do m=1,LIS_rc%nensem(n)
                gid=(t-1)*LIS_rc%nensem(n)+m
                if(  (m.eq.1) &
                     .and. &
                     (count.eq.0) &
                     .and. &
                     ((vardata(gid) .lt. HYMAP2_pe_struc(n)%param_min(i)) &
                     .or. &
                     (vardata(gid) .gt. HYMAP2_pe_struc(n)%param_max(i))) ) then
                   count=count+1   
                   write(LIS_logunit,*) '*****************************************************************', '  ', &
                        'WARNING: HYMAP2 default value is out of LIS-OPT/UE bounds '                , '  ', &
                        'for ', vname                                                             , '  ', &
                        'default value: ', vardata(gid)                                           , '  ', &
                        'parameter min: ', HYMAP2_pe_struc(n)%param_min(i)                        , '  ', &
                        'parameter max: ', HYMAP2_pe_struc(n)%param_max(i)                        , '  ', &   
                        '*****************************************************************'
                   
                endif
             enddo
          enddo
       endif ! if(HYMAP2_pe_struc(n)%param_select(i).eq.1) then 
    enddo  ! do i=1,HYMAP2_pe_struc(n)%nparams
   
    !random initialization
    if(LIS_rc%decSpaceInitMode.eq.1) then  !random initialization 
       seed=seed_base-LIS_localPet !seed must be negative number
       call LIS_rand_func(seed,rand)  !initialize random seed with negative number
       
       do i=1,HYMAP2_pe_struc(n)%nparams
          if(HYMAP2_pe_struc(n)%param_select(i).eq.1) then 
             vname=trim(HYMAP2_pe_struc(n)%param_name(i))
             
             call ESMF_StateGet(DEC_State,vname,varField,rc=status)
             call LIS_verify(status)
             
             call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
                  rc=status)
             call LIS_verify(status)
             
             do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
                do m=1,LIS_rc%nensem(n)
                   if (m.eq.1) then
                      !nothing; leave ensemble 1 with default values
                   else
                      call LIS_rand_func(1,rand)
                      vardata((t-1)*LIS_rc%nensem(n)+m) = &
                           HYMAP2_pe_struc(n)%param_min(i) &
                           + rand * ( HYMAP2_pe_struc(n)%param_max(i) - HYMAP2_pe_struc(n)%param_min(i) )
                   endif
                enddo
             enddo             
          endif
       enddo
    endif

    write(LIS_logunit,*) '[INFO] Finished setting up HYMAP2 decision space '
  end subroutine HYMAP2_setup_pedecvars

end module HYMAP2_peMod
