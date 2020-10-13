!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_AmazonSFobs
! \label{read_AmazonSFobs}
!
! !REVISION HISTORY:
!  2 Oct 2020  Sujay Kumar;   Initial Specification
!
! !INTERFACE: 
subroutine read_AmazonSFobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_fileIOMod
  use AmazonSF_obsMod
  use map_utils

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_Space] Objective Space
!  \end{description}
!
!EOP
  integer                  :: n
  real                     :: time
  integer                  :: offset
  real,    pointer         :: sf(:)
  type(ESMF_Field)         :: sfField
  character*100            :: obsdir
  integer                  :: siteid,grid_index
  type(ESMF_Time)          :: ctime
  logical                  :: data_update
  real                     :: sf2d(LIS_rc%lnc(1),LIS_rc%lnr(1))
  integer                  :: c,r
  integer                  :: status

  n = 1
  sf2d = LIS_rc%udef

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  time = LIS_rc%hr*3600+LIS_rc%mn*60+LIS_rc%ss
  if((mod(time,86400.0).eq.0.0).or.&
       LIS_rc%da.ne.AmazonSF_obs_struc(n)%da) then 

     call ESMF_TimeSet(ctime,yy=LIS_rc%yr,&
          mm=LIS_rc%mo,&
          dd=LIS_rc%da,&
          h=0,&
          m=0,&
          calendar=LIS_calendar,&
          rc=status)
     call LIS_verify(status,'Error in ESMF_TimeSet: readAmazonSFobs')
     AmazonSF_obs_struc(n)%da = LIS_rc%da
     
     offset = nint((ctime-AmazonSF_obs_struc(n)%startTime)/&
          AmazonSF_obs_struc(n)%timestep)+1
     if(offset.gt.0) then 
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(AmazonSF_obs_struc(n)%sites(c,r).gt.0) then 
                 siteid = nint(AmazonSF_obs_struc(n)%sites(c,r))
                 sf2d(c,r) = AmazonSF_obs_struc(n)%sf(siteid,offset)
              endif
           enddo
        enddo
     endif
  endif

  call ESMF_StateGet(Obj_Space,"Amazon_SF",sfField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(sfField,localDE=0,farrayPtr=sf,rc=status)
  call LIS_verify(status)

  sf = LIS_rc%udef
  
  if((mod(time,86400.0).eq.0.0)) then 
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              grid_index =LIS_domain(n)%gindex(c,r)              
              sf(grid_index) = & 
                   sf2d(c,r)
           endif
        enddo
     enddo
  endif

  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)

end subroutine read_AmazonSFobs

!BOP
!
! !ROUTINE: create_AmazonSFobs_filename
! \label(create_AmazonSFobs_filename)
!
! !INTERFACE:
subroutine create_AmazonSFobs_filename(&
     obsdir, &
     yr, mo, nt, fname)
!
! !USES:
  implicit none
!
! !INPUT PARAMETERS:
  character(len=*), intent(in) :: obsdir
  integer  ,        intent(in) :: yr
  integer  ,        intent(in) :: mo
  integer  ,        intent(out) :: nt
  character(len=*), intent(out) :: fname
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
  character (len=4) :: fyr
  integer           :: tmpyr

  if (mo.ge.10) then
     tmpyr = yr + 1
  else
     tmpyr = yr
  endif
  
  write(fyr, '(i4.4)' ) tmpyr
  
  !leap year
  if (((mod(tmpyr,4).eq.0).and.(mod(tmpyr,100).ne.0)).or. &
       (mod(tmpyr,400).eq.0)) then
     nt = 366
  else
     nt = 365
  endif
  
  fname = trim(obsdir)//'/4km_SWE_Depth_WY'//trim(fyr)//'_v01.nc'

end subroutine create_AmazonSFobs_filename


