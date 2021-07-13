!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readNOAAGLDASObs
! \label{readNOAAGLDASObs}
!
! !INTERFACE: 
subroutine readNOAAGLDASObs(source)
! 
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif  
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use NOAAGLDAS_obsMod, only : NOAAGLDASobs

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This subroutine provides the data reader for the
! NOAA GLDAS data  
! 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 July 2021: Sujay Kumar, Initial Specification
! 
!EOP
  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j
  integer           :: ios,nid,tindex,smid,status
  character*100     :: fname
  real              :: sm(NOAAGLDASobs(source)%nc,NOAAGLDASobs(source)%nr,1)
  real              :: sm_data(NOAAGLDASobs(source)%nc*NOAAGLDASobs(source)%nr)
  logical*1         :: sm_data_b(NOAAGLDASobs(source)%nc*NOAAGLDASobs(source)%nr)
  real              :: smobs(LVT_rc%lnc*LVT_rc%lnr)
  logical*1         :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)
  real              :: lat,lon
  type(ESMF_Time)   :: ctime

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  NOAAGLDASobs(source)%smobs = LVT_rc%udef
  smobs= LVT_rc%udef
        
  if(NOAAGLDASobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     NOAAGLDASobs(source)%startmode = .false. 

     fname = NOAAGLDASobs(source)%odir
     inquire(file=trim(fname),exist=file_exists)

     if(file_exists) then
        write(LVT_logunit,*) '[INFO] Reading NOAA GLDAS file ... ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
        ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
        call LVT_verify(ios,'Error opening file '//trim(fname))
        
        ios = nf90_inq_varid(nid, 'soil_moisture_1',smid)
        call LVT_verify(ios, 'Error nf90_inq_varid: soil_moisture_1')

        call ESMF_TimeSet(ctime, yy=LVT_rc%dyr(source), &
          mm = LVT_rc%dmo(source), &
          dd = LVT_rc%dda(source), &
          h = 0, &
          m = 0, &
          calendar = LVT_calendar, &
          rc=status)
             
        tindex = nint((ctime - NOAAGLDASobs(source)%startTime)/&
             NOAAGLDASobs(source)%dt)+1
        if(tindex.gt.0) then

           ios = nf90_get_var(nid, smid, sm,&
                start=(/1,1,tindex/), &
                count=(/NOAAGLDASobs(source)%nc,NOAAGLDASobs(source)%nr,1/))
           call LVT_verify(ios, 'Error nf90_get_var: sm')
           
           ios = nf90_close(ncid=nid)
           call LVT_verify(ios,'Error closing file '//trim(fname))
           
           do r=1, NOAAGLDASobs(source)%nr
              do c=1, NOAAGLDASobs(source)%nc
                 if(sm(c,r,1).lt.0.or.sm(c,r,1).gt.1) then 
                    sm(c,r,1) = LVT_rc%udef
                 endif
              enddo
           enddo
           
           
           do r=1, NOAAGLDASobs(source)%nr
              do c=1, NOAAGLDASobs(source)%nc
                 sm_data(c+(r-1)*NOAAGLDASobs(source)%nc) = sm(c,r,1)
                 if(sm(c,r,1).ne.LVT_rc%udef) then 
                    sm_data_b(c+(r-1)*NOAAGLDASobs(source)%nc) = .true. 
                 else
                    sm_data_b(c+(r-1)*NOAAGLDASobs(source)%nc) = .false.
                 endif
                 if(sm(c,r,1).gt.0.5) then 
                    sm_data(c+(r-1)*NOAAGLDASobs(source)%nc) = LVT_rc%udef
                    sm_data_b(c+(r-1)*NOAAGLDASobs(source)%nc) = .false.
                 endif
              enddo
           enddo
        
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
           call bilinear_interp(LVT_rc%gridDesc(:),&
                sm_data_b, sm_data, smobs_b_ip, smobs, &
                NOAAGLDASobs(source)%nc*NOAAGLDASobs(source)%nr, &
                LVT_rc%lnc*LVT_rc%lnr, &
                NOAAGLDASobs(source)%rlat, NOAAGLDASobs(source)%rlon, &
                NOAAGLDASobs(source)%w11, NOAAGLDASobs(source)%w12, &
                NOAAGLDASobs(source)%w21, NOAAGLDASobs(source)%w22, &
                NOAAGLDASobs(source)%n11, NOAAGLDASobs(source)%n12, &
                NOAAGLDASobs(source)%n21, NOAAGLDASobs(source)%n22, &
                LVT_rc%udef, ios)

        endif
#endif
  
     else
        write(LVT_logunit,*) '[WARN] Unable to read NOAA GLDAS file ...',trim(fname)
        write(LVT_logunit,*) '[WARN]  Note: Check length of filepath, as a long '
        write(LVT_logunit,*) '[WARN] pathname can lead to not reading the file.' 
     endif
     
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smobs(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
              NOAAGLDASobs(source)%smobs(c,r) = smobs(c+(r-1)*LVT_rc%lnc)
           endif
        enddo
     enddo
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       NOAAGLDASobs(source)%smobs,vlevel=1,units="m3/m3")

end subroutine readNOAAGLDASObs

