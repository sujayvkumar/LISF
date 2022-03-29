!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readsyntheticsmObs
! \label{readsyntheticsmObs}
! 
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readsyntheticsmObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod,      only : LDT_rc
  use LDT_timeMgrMod,   only : LDT_get_julss
  use LDT_logMod
  use LDT_DAobsDataMod
  use syntheticsm_obsMod, only : syntheticsmobs
  use map_utils
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the synthetic
! soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  integer           :: ftn
  integer           :: ios,smid
  character*100     :: fname
  real              :: smobs(LDT_rc%lnc(n),LDT_rc%lnr(n))

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LDT_rc%hr)*3600 + 60*LDT_rc%mn + LDT_rc%ss
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  syntheticsmobs(n)%smobs = LDT_rc%udef
  smobs= LDT_rc%udef
        
  if(syntheticsmobs(n)%startmode.or.alarmCheck) then 
     
     syntheticsmobs(n)%startmode = .false. 

     call create_syntheticsm_filename(syntheticsmobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, fname)

     inquire(file=trim(fname),exist=file_exists)
     if(file_exists) then
        write(LDT_logunit,*) 'Reading .. ',trim(fname)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=ftn)
        call LDT_verify(ios,'Error opening file '//trim(fname))
                
        ios = nf90_inq_varid(ftn, 'SoilMoist_tavg',smid)
        call LDT_verify(ios, 'Error nf90_inq_varid: sm')
                
        ios = nf90_get_var(ftn, smid, smobs)
        call LDT_verify(ios, 'Error nf90_get_var: sm')
        
        ios = nf90_close(ncid=ftn)
        call LDT_verify(ios,'Error closing file '//trim(fname))
#endif        


!        ftn = LDT_getNextUnitNumber()
!        open(ftn,file=trim(fname),form='unformatted')
!        read(ftn) smobs
!        call LDT_releaseUnitNumber(ftn)
     endif

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(smobs(c,r).ne.-9999.0) then 
              syntheticsmobs(n)%smobs(c,r) = smobs(c,r)
           endif
        enddo
     enddo
  endif

  call LDT_logSingleDAobs(n,LDT_DAobsData(n)%soilmoist_obs,&
       syntheticsmobs(n)%smobs,vlevel=1)

end subroutine readsyntheticsmObs

!BOP
! !ROUTINE: create_syntheticsm_filename
! \label{create_syntheticsm_filename}
! 
! !INTERFACE: 
subroutine create_syntheticsm_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the synthetic filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the synthetic soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated synthetic filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/SimObs_'//&
       trim(fyr)//trim(fmo)//trim(fda)//'0000.nc'
  
end subroutine create_syntheticsm_filename
