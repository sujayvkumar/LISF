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
! !ROUTINE: get_CRMMScrb
!  \label{get_CRMMScrb}
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
!
! !INTERFACE:
subroutine get_CRMMScrb(n,findex)
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use LIS_timeMgrMod,     only : LIS_tick, LIS_get_nstep
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use CRMMScrb_forcingMod,    only : CRMMScrb_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, CRMMS forcing. 
!  The temporal frequency of the CRMMS data is 3 hours. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer :: t,f
  integer :: ferror
  integer :: tindex
  integer, parameter :: ndays = 10  ! # of days to look back for forcing data
  integer :: order     ! 1 indicates lesser interpolation boundary time
                       ! 2 indicates greater interpolation boundary time
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy2, yr2, mo2, da2, hr2, mn2, ss2
  integer :: movetime  ! 1=move time2 into time1
  real*8  :: timenow, time1, time2
  real*8  :: dumbtime1, dumbtime2
  real    :: gmt1, gmt2,ts1,ts2
  character(len=LIS_CONST_PATH_LEN) :: name
  integer :: nstep

  CRMMScrb_struc(n)%findtime1=0
  CRMMScrb_struc(n)%findtime2=0
  movetime=0
  nstep = LIS_get_nstep(LIS_rc,n)

!-----------------------------------------------------------------
! Determine required CRMMSCRB data times 
! (previous assimilation, current & future assimilation hours)
! The adjustment of the hour and the direction will be done
! in the subroutines that generate the names
!-----------------------------------------------------------------
  yr1=LIS_rc%yr    !Time now
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0        
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr1=LIS_rc%yr    !Previous Hour
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=3*((LIS_rc%hr)/3)
  mn1=0
  ss1=0
  ts1=0
  call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
  
  yr2=LIS_rc%yr    !Next Hour
  mo2=LIS_rc%mo
  da2=LIS_rc%da
  hr2=3*((LIS_rc%hr)/3)
  mn2=0
  ss2=0
  ts2=3*60*60
  call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
  if(timenow.ge.CRMMScrb_struc(n)%CRMMScrbtime2) then
     movetime=1
     CRMMScrb_struc(n)%findtime2=1
  endif
  
  if ( nstep.eq.0 .or. nstep.eq.1 .or.LIS_rc%rstflag(n).eq.1 ) then 
     CRMMScrb_struc(n)%findtime1=1
     CRMMScrb_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
     
!-------------------------------------------------------------------
! Establish CRMMScrbtime1
!-------------------------------------------------------------------
  if (CRMMScrb_struc(n)%findtime1==1) then  
     order=1   
     ferror = 0
     ts1 = -24*60*60
     call CRMMScrbfile(name,CRMMScrb_struc(n)%CRMMScrbdir,&
          yr1,mo1,da1)
     write(LIS_logunit,*) '[INFO] Reading ',trim(name)

     call mapCRMMStimeindex(order, hr1,tindex)

     call read_CRMMScrb(n,name, findex, order,tindex,ferror)
     if ( ferror == 1 ) then 
!-------------------------------------------------------------------
! successfully retrieved forcing data
!-------------------------------------------------------------------
        CRMMScrb_struc(n)%CRMMScrbtime1=time1
     else  
        call LIS_endrun
     end if
  endif
  if(movetime.eq.1) then
     CRMMScrb_struc(n)%CRMMScrbtime1=CRMMScrb_struc(n)%CRMMScrbtime2
     CRMMScrb_struc(n)%findtime2=1 
     do f=1,LIS_rc%met_nf(findex)
        do t=1,LIS_rc%ngrid(n)
           CRMMScrb_struc(n)%metdata1(f,t)=CRMMScrb_struc(n)%metdata2(f,t)
        enddo
     enddo
  endif
  if(CRMMScrb_struc(n)%findtime2.eq.1) then 
     order=2   
     ferror = 0
     ts2 = -24*60*60

     call CRMMScrbfile(name,CRMMScrb_struc(n)%CRMMScrbdir,&
          yr2,mo2,da2)
     write(LIS_logunit,*) '[INFO] Reading ',trim(name)
     call mapCRMMStimeindex(order, hr2,tindex)

     call read_CRMMScrb(n,name,findex,order,tindex,ferror)
     if ( ferror == 1 ) then 
        !-------------------------------------------------------------------
        ! successfully retrieved forcing data
        !-------------------------------------------------------------------
        CRMMScrb_struc(n)%CRMMScrbtime2=time2
     else  
        call LIS_endrun
     end if
  endif
  
end subroutine get_CRMMScrb


!BOP
! !ROUTINE: CRMMScrbfile
! \label{CRMMScrbfile}
!
! !INTERFACE:
subroutine CRMMScrbfile( name, CRMMScrbdir, yr, mo,da)

  implicit none

! !DESCRIPTION:
!   This subroutine puts together a filename for
!   a given timestamp
! 
!  The arguments are:
!  \begin{description}
!  \item[CRMMScrbdir]
!    Name of the CRMMS directory
!  \item[yr]
!    year 
!  \item[mo]
!    month of year
!  \item[da]
!    month of year
!   \item[name]
!   name of the CRMMScrb  file
!  \end{description}
!
!EOP

! !ARGUMENTS: 
  integer :: yr, mo,da

  character(len=*) :: name
  character(len=*) :: CRMMScrbdir
  character(4) :: cyear
  character(2) :: cmo,cda

   write ( cyear, '(i4)' ) yr
   write ( cmo, '(i2.2)' ) mo
   write ( cda, '(i2.2)' ) da

   name = trim(CRMMScrbdir)//'/mm_'//trim(cyear)//&
        trim(cmo)//trim(cda)//'.nc'

end subroutine CRMMScrbfile

subroutine mapCRMMStimeindex(order, hr,tindex)

  implicit none

  integer     :: order
  integer     :: hr
  integer     :: tindex


  if(hr.ge.0.and.hr.le.2) then
     tindex = 1
  elseif(hr.ge.3.and.hr.le.5) then
     tindex = 2
  elseif(hr.ge.6.and.hr.le.8) then
     tindex = 3
  elseif(hr.ge.9.and.hr.le.11) then
     tindex = 4
  elseif(hr.ge.12.and.hr.le.14) then
     tindex = 5
  elseif(hr.ge.15.and.hr.le.17) then
     tindex = 6
  elseif(hr.ge.18.and.hr.le.20) then
     tindex = 7
  elseif(hr.ge.21.and.hr.le.23) then
     tindex = 8        
  endif

end subroutine mapCRMMStimeindex
