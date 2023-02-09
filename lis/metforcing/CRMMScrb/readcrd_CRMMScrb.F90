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
! !ROUTINE: readcrd_CRMMScrb
!  \label{readcrd_CRMMScrb}
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
!
! !INTERFACE:    
subroutine readcrd_CRMMScrb()
! !USES:
  use ESMF
  use LIS_coreMod,       only : LIS_rc, LIS_config
  use LIS_logMod,        only : LIS_logunit, LIS_verify
  use CRMMScrb_forcingMod, only : CRMMScrb_struc

  implicit none

! !DESCRIPTION:
!
!  This routine reads the options specific to CRMMS forcing data
!  
!EOP

  integer :: n,rc
  
  call ESMF_ConfigFindLabel(LIS_config,"CRMMS met forcing directory:",rc=rc)
  call LIS_verify(rc, 'CRMMS met forcing directory: not defined')
  do n=1,LIS_rc%nnest    
     call ESMF_ConfigGetAttribute(LIS_config,CRMMScrb_struc(n)%CRMMScrbdir,rc=rc)
  enddo


  write(unit=LIS_logunit,fmt=*)'[INFO] Using CRMMS met forcing'

  do n=1,LIS_rc%nnest
     write(unit=LIS_logunit,fmt=*) '[INFO] CRMMS met forcing directory :',trim(CRMMScrb_struc(n)%CRMMScrbdir)

     CRMMScrb_struc(n)%CRMMScrbtime1 = 3000.0
     CRMMScrb_struc(n)%CRMMScrbtime2 = 0.0
  enddo

end subroutine readcrd_CRMMScrb
