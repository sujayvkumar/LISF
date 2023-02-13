!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_CRMMScrb
!  \label{read_CRMMScrb}
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
! 
! !INTERFACE:
subroutine read_CRMMScrb(n, fname, findex, order, tindex,ferror)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod, only : LIS_forc
  use LIS_logMod
  use CRMMScrb_forcingMod, only : CRMMScrb_struc
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: findex
  integer, intent(in)          :: order
  character(len=*), intent(in) :: fname
  integer                      :: tindex
  integer, intent(out)         :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  output files, transforms into 2 LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!EOP
  integer                 :: c,r,gindex
  integer                 :: ftn
  integer                 :: vindex
  logical                 :: file_exists
  integer                 :: t2Id, pcpId
  real                    :: t2in(CRMMScrb_struc(n)%ncol,CRMMScrb_struc(n)%nrow)
  real                    :: pcpin(CRMMScrb_struc(n)%ncol,CRMMScrb_struc(n)%nrow)

  real                    :: t2out(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                    :: pcpout(LIS_rc%lnc(n),LIS_rc%lnr(n))  

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  inquire (file=trim(fname), exist=file_exists)
  if (file_exists) then      
     ferror = 1
     call LIS_verify(nf90_open(path=trim(fname),mode=NF90_NOWRITE, &
          ncid = ftn), 'nf90_open failed in read_CRMMScrb')
     call LIS_verify(nf90_inq_varid(ftn,'otgrid',t2Id),&
          'nf90_inq_varid failed for temperature in read_CRMMScrb')
     call LIS_verify(nf90_inq_varid(ftn,'ogrid',pcpId),&
          'nf90_inq_varid failed for precip in read_CRMMScrb')
     
     call LIS_verify(nf90_get_var(ftn,t2id,t2in,&
          start=(/1,1,tindex/),&
          count=(/CRMMScrb_struc(n)%ncol,CRMMScrb_struc(n)%nrow,1/)),&
          'nf90_get_var failed for t2 in read_CRMMScrb')
     call LIS_verify(nf90_get_var(ftn,pcpid,pcpin,&
          start=(/1,1,tindex/),&
          count=(/CRMMScrb_struc(n)%ncol,CRMMScrb_struc(n)%nrow,1/)),&
          'nf90_get_var failed for precip in read_CRMMScrb')
     
     call interp_CRMMScrb( n, findex, 1, t2in, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n), LIS_rc%lnr(n), t2out )
     call interp_CRMMScrb( n, findex, 2, pcpin, LIS_rc%gridDesc(n,:), &
          LIS_rc%lnc(n), LIS_rc%lnr(n), pcpout )     

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 CRMMScrb_struc(n)%metdata1(1,&
                      LIS_domain(n)%gindex(c,r)) =&
                      t2out(c,r)
                 CRMMScrb_struc(n)%metdata1(2,&
                      LIS_domain(n)%gindex(c,r)) =&
                      pcpout(c,r)                 
              else
                 CRMMScrb_struc(n)%metdata2(1,&
                      LIS_domain(n)%gindex(c,r)) = &
                      t2out(c,r)
                 CRMMScrb_struc(n)%metdata2(2,&
                      LIS_domain(n)%gindex(c,r)) =&
                      pcpout(c,r)                 
                 
              endif
           endif
        enddo
     enddo
           
  else
     write(LIS_logunit,*) '[ERR] Forcing file '//trim(fname)//' not found'
     ferror = 0 
  endif
#endif

  
end subroutine read_CRMMScrb

!BOP
!
! !ROUTINE: interp_CRMMScrb
! \label{interp_CRMMScrb}
!
! !INTERFACE:
subroutine interp_CRMMScrb (n, findex,  vindex, var_in, lis_gds, &
                        nc, nr, varfield )

! !USES:
  use CRMMScrb_forcingMod, only : CRMMScrb_struc
  use LIS_coreMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: vindex
  integer                :: nc      ! Number of columns (in the E-W dimension) in the LIS grid
  integer                :: nr      ! Number of rows (in the N-S dimension) in the LIS grid
  real                   :: var_in(CRMMScrb_struc(n)%ncol,CRMMScrb_struc(n)%nrow)
  real                   :: lis_gds(50)
  real                   :: varfield(nc,nr)

! !DESCRIPTION:
!   This subroutine interpolates a given CRMMS field
!   to the LIS grid.
!  
!EOP

  integer :: iret
  integer :: mo
  logical*1              :: lb(CRMMScrb_struc(n)%ncol*CRMMScrb_struc(n)%nrow)
  integer :: count1, i, j,c,r

  real    :: var1d(CRMMScrb_struc(n)%ncol*CRMMScrb_struc(n)%nrow)
  real    :: lis1d(nc*nr)   ! LIS Grid
  logical*1 :: lo(nc*nr)    ! LIS Grid

!=== End variable declarations

!--------------------------------------------------------------------
! NOTE:: Recommended to use budget bilinear for precip forcing fields
!--------------------------------------------------------------------
  CRMMScrb_struc(n)%mi = CRMMScrb_struc(n)%ncol*CRMMScrb_struc(n)%nrow
  mo = nc * nr   ! LIS Grid

  lo = .true.    !-Initialize output bitmap
  lb = .false.  


  do r=1,CRMMScrb_struc(n)%nrow
     do c=1,CRMMScrb_struc(n)%ncol
        if(CRMMScrb_struc(n)%mask(c,r).eq.1) then 
           var1d(c+(r-1)*CRMMScrb_struc(n)%ncol) = &
                var_in(c,r)
        else
           var1d(c+(r-1)*CRMMScrb_struc(n)%ncol) = -9999.0
        endif
        if(var1d(c+(r-1)*CRMMScrb_struc(n)%ncol).ne.-9999.0) then 
           lb(c+(r-1)*CRMMScrb_struc(n)%ncol) = .true.
        endif
     enddo
  enddo
  
  if ( LIS_isatAfinerResolution(n,0.008333) ) then
     if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
        ! == BILINEAR INTERPOLATION ====
        call bilinear_interp( lis_gds,  lb, var1d, lo, &
             lis1d, CRMMScrb_struc(n)%mi, mo, &
             LIS_domain(n)%lat, LIS_domain(n)%lon,&
             CRMMScrb_struc(n)%w111, CRMMScrb_struc(n)%w121, &
             CRMMScrb_struc(n)%w211, CRMMScrb_struc(n)%w221, &
             CRMMScrb_struc(n)%n111, CRMMScrb_struc(n)%n121, &
             CRMMScrb_struc(n)%n211, CRMMScrb_struc(n)%n221, LIS_rc%udef, iret)
     endif
     
  else
     call upscaleByAveraging(CRMMScrb_struc(n)%mi, mo, LIS_rc%udef, &
          CRMMScrb_struc(n)%n111, lb, var_in, lo, lis1d)
  endif
!--------------------------------------------------------------------
! Create 2D array (on LIS LIS_domain) for main program.
!--------------------------------------------------------------------
  if(vindex.eq.1) then 
     count1 = 0
     do j = 1, nr
        do i = 1, nc
           if(lo(i+count1)) then 
              varfield(i,j) = (lis1d(i+count1)-32)*5/9+273.15
           else
              varfield(i,j) = -9999.0
           endif
        enddo
        count1 = count1 + nc
     enddo
  elseif(vindex.eq.2) then
     count1 = 0
     do j = 1, nr
        do i = 1, nc
           if(lo(i+count1)) then 
              varfield(i,j) = lis1d(i+count1)*25.4  !inches to mm
           else
              varfield(i,j) = -9999.0
           endif
        enddo
        count1 = count1 + nc
     enddo
  endif

   
end subroutine interp_CRMMScrb
