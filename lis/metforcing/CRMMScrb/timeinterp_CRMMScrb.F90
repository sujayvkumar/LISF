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
! !ROUTINE: timeinterp_CRMMScrb
! \label{timeinterp_CRMMScrb}
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
!
! !INTERFACE:
subroutine timeinterp_CRMMScrb(n, findex)
! !USES:
  use ESMF
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc, LIS_FORC_Base_State
  use LIS_coreMod,       only : LIS_rc,LIS_domain
  use LIS_constantsMod,  only : LIS_CONST_SOLAR
  use LIS_timeMgrMod,    only : LIS_tick, LIS_time2date
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use CRMMScrb_forcingMod, only : CRMMScrb_struc
 
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing source
!  \end{description}
! 
!EOP
  integer          :: zdoy
  real             :: zw1, zw2
  real             :: czm, cze, czb
  real             :: wt1, wt2,swt1,swt2
  real             :: gmt1, gmt2
  integer          :: t,index1
  integer          :: bdoy,byr,bmo,bda,bhr,bmn
  real*8           :: btime,newtime1,newtime2
  real             :: tempgmt1,tempgmt2
  integer          :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer          :: tempbss
  real             :: tempbts
  integer          :: status
  type(ESMF_Field) :: tmpField
  type(ESMF_Field) :: pcpField
  real,pointer     :: tmp(:)
  real,pointer     :: pcp(:)
  logical          :: forcing_z, forcing_ch
  
  btime=CRMMScrb_struc(n)%CRMMScrbtime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)
  
  tempbdoy=bdoy
  tempgmt1=gmt1
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,& 
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, & 
       tempbss,tempbts)
  
  btime=CRMMScrb_struc(n)%CRMMScrbtime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy=bdoy
  tempgmt2=gmt2
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)
  
!=== Interpolate Data in time      
  wt1=(CRMMScrb_struc(n)%CRMMScrbtime2-LIS_rc%time)/ & 
      (CRMMScrb_struc(n)%CRMMScrbtime2-CRMMScrb_struc(n)%CRMMScrbtime1)
  wt2=1.0-wt1
  swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Tair in the forcing variables list')


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable Rainf in the forcing variables list')


  !do block precipitation interpolation
  !WRF writes precipitation accumlated from the beginning of the
  !simulation.  Therefore compute the difference between the bookends
  !to find the accumulated precipitation for the hour.
  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if (CRMMScrb_struc(n)%metdata2(2,index1).ne.LIS_rc%udef) then 
        pcp(t) = CRMMScrb_struc(n)%metdata2(2,index1) 
        pcp(t) = pcp(t)/(3*60.0*60.0)
     endif
  enddo

  !linearly interpolate everything else

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index 
     if((CRMMScrb_struc(n)%metdata1(1,index1).ne.LIS_rc%udef) .and. &
        (CRMMScrb_struc(n)%metdata2(1,index1).ne.LIS_rc%udef)) then 
        tmp(t) = CRMMScrb_struc(n)%metdata1(1,index1)*wt1 + & 
                 CRMMScrb_struc(n)%metdata2(1,index1)*wt2
     endif
  enddo

end subroutine timeinterp_CRMMScrb
 
