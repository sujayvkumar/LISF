!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_set_pedecvars
!  \label{HYMAP2_set_pedecvars}
!
! !REVISION HISTORY:
! 02 Oct 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine HYMAP2_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use HYMAP2_routingMod, only : HYMAP2_routing_struc
  use HYMAP2_peMod,  only : HYMAP2_pe_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to HYMAP2 model variables. 
! 
!EOP
  integer                :: n
  real, pointer          :: vdata(:)
  character*100          :: vname
  integer, pointer       :: mod_flag_HYMAP2(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_HYMAP2(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_HYMAP2 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,HYMAP2_pe_struc(n)%nparams
     if(HYMAP2_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(HYMAP2_pe_struc(n)%param_name(i))
        call HYMAP2_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call HYMAP2_checkBounds(n,DEC_State,vname, vdata, mod_flag_HYMAP2)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call HYMAP2_checkConstraints(n,DEC_State, mod_flag_HYMAP2)

  !set variables given modflag; if flag set will leave values alone
  call HYMAP2_setVars(n,DEC_State,mod_flag_HYMAP2)

  !send mod flag to ESMF state (feasibility flag)
  call HYMAP2_setModFlag(n,DEC_State,Feas_State,mod_flag_HYMAP2)
end subroutine HYMAP2_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine HYMAP2_getvardata(n,DEC_State,vname, vdata, statusStateGet)
! !USES:
  use ESMF
  use HYMAP2_routingMod
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real          :: vdata(HYMAP2_routing_struc(n)%nseqall)
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
! 
!EOP
  real, pointer          :: vardata(:)
  type(ESMF_Field)       :: varField
  integer                :: statusStateGet, statusFieldGet,i
  
  call ESMF_StateGet(DEC_State,vname,varField,rc=statusStateGet)
!  call LIS_verify(status)
  
  if(statusStateGet.eq.0) then
     call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
          rc=statusFieldGet)
     call LIS_verify(statusFieldGet)
     vdata=vardata
  endif
  
end subroutine HYMAP2_getvardata

subroutine HYMAP2_checkBounds(n,DEC_State,vname, vardata, mod_flag_HYMAP2)
! !USES:
  use ESMF
  use HYMAP2_routingMod
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real          :: vardata(HYMAP2_routing_struc(n)%nseqall)
  integer       :: mod_flag_HYMAP2(HYMAP2_routing_struc(n)%nseqall)
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  integer                :: status
  integer                :: t

  call ESMF_StateGet(DEC_State,vname,varField,rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
  call LIS_verify(status)
  
  do t=1,HYMAP2_routing_struc(n)%nseqall
     if(vardata(t).lt.vardata_min) then 
        mod_flag_HYMAP2(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_HYMAP2(t) = 1
     endif
  enddo
end subroutine HYMAP2_checkBounds

subroutine HYMAP2_checkConstraints(n,DEC_State,mod_flag_HYMAP2)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use HYMAP2_routingMod, only : HYMAP2_routing_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_HYMAP2(HYMAP2_routing_struc(n)%nseqall)
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to HYMAP2 model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  character*100          :: vname
  integer                :: t
  integer       :: status1, status2
  real, allocatable :: vardata1(:)
  real, allocatable :: vardata2(:)

  !nothing to check

end subroutine HYMAP2_checkConstraints

subroutine HYMAP2_setVars(n,DEC_State,mod_flag_HYMAP2)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use HYMAP2_routingMod, only : HYMAP2_routing_struc
  use HYMAP2_peMod,  only : HYMAP2_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_HYMAP2(HYMAP2_routing_struc(n)%nseqall)
  type(ESMF_State)       :: DEC_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real          :: vardata(HYMAP2_routing_struc(n)%nseqall)
  character*100          :: vname
  integer                :: i,t,k,m, status

  do i=1,HYMAP2_pe_struc(n)%nparams
     if(HYMAP2_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(HYMAP2_pe_struc(n)%param_name(i))
        call HYMAP2_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,HYMAP2_routing_struc(n)%nseqall
           do m=1,LIS_rc%nensem(n)
              k=(t-1)*LIS_rc%nensem(n)+m
              if(mod_flag_HYMAP2(k).eq.0) then 
                 if(vname.eq."RIVWTH") &
                      HYMAP2_routing_struc(n)%rivwth(t,m) = vardata(k) 
                 if(vname.eq."RIVHGT") & 
                      HYMAP2_routing_struc(n)%rivhgt(t,m) = vardata(k) 
                 if(vname.eq."RIVROUGHNESS") &
                      HYMAP2_routing_struc(n)%rivman(t,m) = vardata(k) 
              endif
           enddo
        enddo
     endif
  enddo
end subroutine HYMAP2_setVars

subroutine HYMAP2_setModFlag(n,DEC_State,Feas_State,mod_flag_HYMAP2)
! !USES:
  use ESMF
  use HYMAP2_routingMod
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer                :: mod_flag_HYMAP2(HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n))
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine sets the feasibility flag
! 
!EOP
  type(ESMF_Field)       :: feasField
  integer                :: t
  integer                :: status
  integer, pointer       :: modflag(:)

  call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(feasField,localDE=0,farrayPtr=modflag,rc=status)
  call LIS_verify(status)

  do t=1,HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n)
     if(mod_flag_HYMAP2(t).eq.1) then 
        modflag(t)=1
     endif
  enddo

end subroutine HYMAP2_setModFlag
