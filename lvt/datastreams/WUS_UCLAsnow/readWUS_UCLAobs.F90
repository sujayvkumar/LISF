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
! 
! !ROUTINE: readWUS_UCLAObs
! \label{readWUS_UCLAObs}
!
! !INTERFACE:
subroutine readWUS_UCLAObs(source)
! 
  ! !USES:
  use netcdf
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use WUS_UCLA_obsMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for WUS_UCLA SWE and snowdepth data.
! LVT expects the data to be provided in a timestamped manner. The raw
! WUS_UCLA binary files are read, and the data is interpolated to the
! LIS model output. The data is interpolated using the bilinear
! averaging technique.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 May 2010: Sujay Kumar, Initial Specification
!  1 Aug 2019: David Mocko, Added/fixed log messages
! 
!EOP
  integer                :: source
  integer                :: yr, mo, da, hr
  integer                :: i,j,t,c,r
  integer                :: stn_col, stn_row
  real                   :: col,row
  character*100          :: fname
  logical                :: file_exists
  logical                :: readflag
  integer                :: ftn, ios,sndid
  integer                :: status
  type(ESMF_Time)        :: snodastime, snodastime1
  integer                :: stnindex,tind
  real                   :: offset
  logical                :: alarmCheck
  integer                :: iret
  real                   :: timenow
  real                   :: udef
  logical*1              :: lb(wusuclaobs(source)%nc*wusuclaobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: var(wusuclaobs(source)%nc, wusuclaobs(source)%nr)
  real                   :: var1d(wusuclaobs(source)%nc*wusuclaobs(source)%nr)
  real                   :: var_ip(LVT_rc%lnc*LVT_rc%lnr)
  
!process snow depth data
  var_ip  = LVT_rc%udef
  udef = -9999.0

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
     LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)
  if(wusuclaobs(source)%startflag.or.alarmCheck) then

     wusuclaobs(source)%startflag  = .false. 

     call create_WUS_UCLA_snwdfilename(wusuclaobs(source)%odir, &
          LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), &
          fname)
     
     inquire(file=trim(fname),exist=file_exists)

     if (file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading WUS_UCLA SNWD data ', &
                             trim(fname)
        call LVT_verify(nf90_open(path=fname,mode=NF90_NOWRITE,ncid=ftn),&
             'Error opening WUS_UCLA file')
        call LVT_verify(nf90_inq_varid(ftn,'SD_Post',sndId),&
             'nf90_inq_varid failed for SD_Post in read_WUS_UCLAobs')
        call LVT_verify(nf90_get_var(ftn,sndId,var,&
             count=(/wusuclaobs(source)%nc,wusuclaobs(source)%nr,1/),&
             start=(/wusuclaobs(source)%c_off,wusuclaobs(source)%r_off,1/)),&
             'nf90_get_var failed for SD_Post in read_WUS_UCLAobs')
        call LVT_verify(nf90_close(ftn),&
             'nf90_close failed in readWUS_UCLAobs')
        
        lb = .false.
        var1d = -9999.0
        do r=1,wusuclaobs(source)%nr
           do c=1,wusuclaobs(source)%nc
              if(var(c,r).ge.0) then 
                 var1d(c+(r-1)*wusuclaobs(source)%nc) = var(c,r)
                 lb(c+(r-1)*wusuclaobs(source)%nc) = .true. 
              endif
           enddo
        enddo

        call upscaleByAveraging(&
             wusuclaobs(source)%nc*wusuclaobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,&
             LVT_rc%udef,&
             wusuclaobs(source)%n11,&
             lb,var1d,lo,var_ip)

        readflag = .false. 

        write(LVT_logunit,*) '[INFO] Successfully processed SNWD from ', &
                             trim(fname)
     else
        write(LVT_logunit,*) '[WARN] Could not find WUS_UCLA SNWD file ', &
                             trim(fname)
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH,source,var_ip,vlevel=1,units="m")

end subroutine readWUS_UCLAObs


!BOP
! 
! !ROUTINE: create_WUS_UCLA_snwdfilename
! \label{create_WUS_UCLA_snwdfilename}
!
! !INTERFACE:
subroutine create_WUS_UCLA_snwdfilename(odir, yr,mo,da, fname)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in)  :: odir
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo
  integer,          intent(in)  :: da
! 
! !OUTPUT PARAMETERS:
  character(len=*), intent(out) :: fname
!
! !DESCRIPTION: 
! This routine creates a timestamped WUS_UCLA filename.
! The filename convention is as follows: 
!   us : region (united states)
!   ssm : simple snow model
!   v1  : operational snow model output
!   1034 : Snow water equivalent
!   ts__ : integral through all layers of the snow pack
!   T0001 : A 1-hour snapshot
! 
!  The arguments are: 
!  \begin{description}
!   \item[odir] WUS_UCLA base directory
!   \item[yr]   year of data 
!   \item[mo]   month of data
!   \item[da]   day of data
!   \item[fname]  Name of the WUS_UCLA file  
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4             :: fyr
  character*2             :: fmo
  character*2             :: fda
  
  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo
  write(fda, '(i2.2)' ) da

  fname = trim(odir)//'/'//trim(fyr)//&
       '/WUS_UCLA_SR_v01_agg_16_SWE_SCA_SD_POST_' &
       //fyr//fmo//fda//'.nc4'
  
end subroutine create_WUS_UCLA_snwdfilename
