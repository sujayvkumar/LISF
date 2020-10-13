!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readAmazonSFObs
! \label{readAmazonSFObs}
!
! !INTERFACE: 
subroutine readAmazonSFObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,      only : LVT_rc, LVT_domain
  use LVT_logMod,       only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use LVT_timeMgrMod,   only : LVT_calendar, LVT_tick
  use AmazonSF_obsMod,  only : amazonsfobs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer,   intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for AmazonSF station data. 
! The plugin processes latent, sensible, and ground heat flux data
! from the in-situ measurements. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Apr 2020: Sujay Kumar, Initial Specification
! 
!EOP
!BOP
! !ARGUMENTS: 

!EOP
  type(ESMF_Time)  :: ctime
  integer          :: t,c,r
  real             :: gmt
  real                :: time
  integer             :: status
  integer             :: offset
  integer             :: siteid
  integer             :: k
  real                :: sf(LVT_rc%lnc,LVT_rc%lnr)

  sf  = LVT_rc%udef 

  time = LVT_rc%dhr(source)*3600+LVT_rc%dmn(source)*60+LVT_rc%dss(source)
  if((mod(time,86400.0).eq.0.0).or.&
       (LVT_rc%dda(source).ne.amazonsfobs(source)%da)) then

     print*, LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source)
     call ESMF_TimeSet(ctime, yy=LVT_rc%dyr(source),&
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
     call LVT_verify(status,'error in timeset: readArmobs(Source)')
     amazonsfobs(source)%da = LVT_rc%dda(source)

     offset = nint((ctime-amazonsfobs(source)%startTime)/&
          amazonsfobs(source)%ts) + 1

     if(offset.gt.0) then 
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(amazonsfobs(source)%sites(c,r).gt.0) then 
                 siteid = nint(amazonsfobs(source)%sites(c,r))
                 sf(c,r) = amazonsfobs(source)%sf(siteid,offset)
!                 print*, LVT_rc%dyr(source),LVT_rc%dmo(source),&
!                      LVT_rc%dda(source),&
!                      c,r,sf(c,r)
              endif
           enddo
        enddo        
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_streamflow,source,sf,&
       vlevel=1,units="m3/s")

  
end subroutine readAmazonSFObs

