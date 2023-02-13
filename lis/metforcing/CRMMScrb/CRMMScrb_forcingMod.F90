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
module CRMMScrb_forcingMod
!BOP
! !MODULE: CRMMScrb_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the forcing data used in the Colorado
!  River Mid-term Modeling System (CRMMS). This data contains air
!  temperature and precipitation inputs
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif
  
  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_CRMMScrb      !defines the native resolution of 
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: CRMMScrb_struc
!EOP

  type, public    :: CRMMScrb_type_dec
     character(len=LIS_CONST_PATH_LEN) :: CRMMScrbdir
     real         :: ts
     real*8       :: CRMMScrbtime1,CRMMScrbtime2
     integer      :: findtime1,findtime2
     integer      :: startflag
     integer                :: ncol,nrow,mi
     integer, allocatable   :: n111(:), n121(:)
     integer, allocatable   :: n211(:), n221(:)
     real, allocatable      :: w111(:), w121(:)
     real, allocatable      :: w211(:), w221(:)
     
     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:)
     character*100     :: maskfile
     real, allocatable :: mask(:,:)

  end type CRMMScrb_type_dec

  type(CRMMScrb_type_dec), allocatable :: CRMMScrb_struc(:)
!EOP
contains
  
!BOP
!
! !ROUTINE: init_CRMMScrb
! \label{init_CRMMScrb}
!
! !REVISION HISTORY: 
! 
! !INTERFACE:
  subroutine init_CRMMScrb(findex)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    implicit none
! !ARGUMENTS:  
    integer, intent(in) :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for CRMMS 
!  forcing data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!  Note that no interpolation is performed for
!  this forcing. The data is expected to be in the same map projection
!  and resolution as that of the current LIS run.
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing source
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_CRMMScrb](\ref{readcrd_CRMMScrb}) \newline
!     reads the runtime options specified for WRF output data
!  \end{description}
!EOP
    
    integer :: n
    integer :: ftn,maskid
    logical :: file_exists
    real                   :: gridDesci(50)
    
    
    allocate(CRMMScrb_struc(LIS_rc%nnest))

    call readcrd_CRMMScrb()

    LIS_rc%met_nf(findex) = 2

    do n=1, LIS_rc%nnest

       allocate(CRMMScrb_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(CRMMScrb_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       CRMMScrb_struc(n)%metdata1 = 0
       CRMMScrb_struc(n)%metdata2 = 0


       CRMMScrb_struc(n)%ts = 3*60*60 
       call LIS_update_timestep(LIS_rc, n, CRMMScrb_struc(n)%ts)

       gridDesci = 0 
       gridDesci(1) = 0
       gridDesci(2) = 1321
       gridDesci(3) = 1681
       gridDesci(4) = 29.9999995
       gridDesci(5) = -116.0000005
       gridDesci(6) = 128
       gridDesci(7) = 43.9994395
       gridDesci(8) = -105.0004405
       gridDesci(9) = 0.00833
       gridDesci(10) = 0.00833
       gridDesci(20) = 64

       CRMMScrb_struc(n)%ncol = 1321
       CRMMScrb_struc(n)%nrow = 1681

       
       CRMMScrb_struc(n)%mi = CRMMScrb_struc(n)%ncol * CRMMScrb_struc(n)%nrow

       CRMMScrb_struc(n)%startFlag = .true. 
       
       if ( LIS_isatAfinerResolution(n,gridDesci(9)) ) then
          if (trim(LIS_rc%met_interp(findex)) .eq. "bilinear") then
             allocate(CRMMScrb_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             allocate(CRMMScrb_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
             
             call bilinear_interp_input(n, gridDesci(:),&
                  CRMMScrb_struc(n)%n111, CRMMScrb_struc(n)%n121, &
                  CRMMScrb_struc(n)%n211, CRMMScrb_struc(n)%n221, &
                  CRMMScrb_struc(n)%w111, CRMMScrb_struc(n)%w121, &
                  CRMMScrb_struc(n)%w211, CRMMScrb_struc(n)%w221 )
          endif
       else
             
          allocate(CRMMScrb_struc(n)%n111(CRMMScrb_struc(n)%mi))

          call upscaleByAveraging_input(&
               gridDesci(:),              &
               LIS_rc%gridDesc(n,:),        &
               CRMMScrb_struc(n)%mi,            &
               LIS_rc%lnc(n)*LIS_rc%lnr(n), &
               CRMMScrb_struc(n)%n111)

       endif

       allocate(CRMMScrb_struc(n)%mask(CRMMScrb_struc(n)%ncol,&
            CRMMScrb_struc(n)%nrow))
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
       inquire (file=trim(CRMMScrb_struc(n)%maskfile), exist=file_exists)
       if (file_exists) then      
          call LIS_verify(nf90_open(path=trim(CRMMScrb_struc(n)%maskfile),&
               mode=NF90_NOWRITE, &
               ncid = ftn), 'nf90_open failed init_CRMMScrb')
          call LIS_verify(nf90_inq_varid(ftn,'CRMMS_mask',maskId),&
               'nf90_inq_varid failed for CRMMS_mask in init_CRMMScrb')
          
          call LIS_verify(nf90_get_var(ftn,maskid,CRMMScrb_struc(n)%mask),&
          'nf90_get_var failed for CRMMS_mask in init_CRMMScrb')
          call LIS_verify(nf90_close(ftn))
          
       else
          CRMMScrb_struc(n)%mask = -9999.0
       endif
#endif       
    enddo

  end subroutine init_CRMMScrb
end module CRMMScrb_forcingMod
