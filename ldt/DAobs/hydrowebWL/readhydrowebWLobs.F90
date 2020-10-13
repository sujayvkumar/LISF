!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readhydrowebWLobs
! \label{readhydrowebWLobs}
! 
! !REVISION HISTORY: 
!  28 May 2019: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readhydrowebWLobs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_DAobsDataMod
  use hydrowebWL_obsMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the 
! LPRM vegetation optical depth retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r
  integer           :: status
  integer           :: ftn
  integer           :: tId,wlId,dim1Id, dim2Id
  type(ESMF_Time)   :: startTime, currentTime
  integer           :: offset
  character*100     :: fname
  real              :: wlobs(LDT_rc%lnc(n),LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------
  wlobs= LDT_rc%udef


  inquire(file=hydrowebWLobs(n)%odir, exist=file_exists) 

  if(file_exists.and.hydrowebWLobs(n)%startFlag) then
     
     hydrowebWLobs(n)%startFlag = .false.
     write(LDT_logunit,*) '[INFO] Reading ',trim(hydrowebWLobs(n)%odir)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     status = nf90_open(path=trim(hydrowebWLobs(n)%odir),&
          mode=NF90_NOWRITE,ncid=ftn)
     call LDT_verify(status,'Error opening file '//trim(hydrowebWLobs(n)%odir))
     
     call LDT_verify(nf90_inq_dimid(ftn,'time',dim1Id),&
          'Error with nf90_inq_dimid: time')
     call LDT_verify(nf90_inq_dimid(ftn,'sites',dim2Id),&
          'Error with nf90_inq_dimid: sites')

     call LDT_verify(nf90_inquire_dimension(ftn,dim1id,&
          len=hydrowebWLobs(n)%ntimes),&
          'Error with nf90_inquire_dimension: time')
     call LDT_verify(nf90_inquire_dimension(ftn,dim2id,&
          len=hydrowebWLobs(n)%nsites),&
          'Error with nf90_inquire_dimension: sites')

     allocate(hydrowebwlobs(n)%wlobs(&
          hydrowebwlobs(n)%nsites,&
          hydrowebwlobs(n)%ntimes))

     allocate(hydrowebwlobs(n)%time(&
          hydrowebwlobs(n)%ntimes))

     call LDT_verify(nf90_inq_varid(ftn,'time',tId),&
          'Error with nf90_inq_varid: time')

     call LDT_verify(nf90_get_var(ftn,tId,hydrowebwlobs(n)%time), &
          'Error with nf90_get_var: time')

     call LDT_verify(nf90_inq_varid(ftn,'water_elevation',wlId),&
          'Error with nf90_inq_varid: water_elevation')
     
     call LDT_verify(nf90_get_var(ftn,wlId, hydrowebwlobs(n)%wlobs),&
          'Error with nf90_get_var: water_elevation')

     status = nf90_close(ncid=ftn)
     call LDT_verify(status,'Error closing file '//trim(fname))
#endif
  endif

  timenow = float(LDT_rc%hr*3600 + 60*LDT_rc%mn +LDT_rc%ss)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  wlobs = -9999.0
  
  if(alarmCheck) then 
     
     call ESMF_TimeSet(startTime, yy=2016,&
          mm = 7, dd=1, h =0, &
          m = 0, s = 0, calendar = LDT_calendar,&
          rc=status)
     call LDT_verify(status, 'ESMF_TimeSet failed in readhydrowebWLobs')
     
     call ESMF_TimeSet(currentTime, yy=LDT_rc%yr,&
          mm = LDT_rc%mo, dd=LDT_rc%da, h =LDT_rc%hr, &
          m = LDT_rc%mn, s = 0, calendar = LDT_calendar,&
          rc=status)
     call LDT_verify(status, 'ESMF_TimeSet failed in readhydrowebWLobs')
     
     offset = nint((currentTime - startTime)/hydrowebwlobs(n)%ts) + 1
     if(offset.gt.0) then 
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(hydrowebwlobs(n)%sites_data(c,r).gt.0) then 
                 wlobs(c,r) = hydrowebwlobs(n)%wlobs(&
                      nint(hydrowebwlobs(n)%sites_data(c,r)),offset)
              endif
           enddo
        enddo
     endif
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%wl_obs,&
       wlobs,vlevel=1)

end subroutine readhydrowebWLobs

