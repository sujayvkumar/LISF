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
! !MODULE: reset_CRMMScrb
!  \label{reset_CRMMScrb}
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
! 
! !INTERFACE:
subroutine reset_CRMMScrb()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use CRMMScrb_forcingMod, only : CRMMScrb_struc
!
! !DESCRIPTION:
!  Routine to reset CRMMScrb forcing related memory allocations.   
! 
!EOP
  implicit none
  
  integer   :: n

  do n=1,LIS_rc%nnest
     CRMMScrb_struc(n)%CRMMScrbtime1 = 3000.0
     CRMMScrb_struc(n)%CRMMScrbtime2 = 0.0
  enddo

end subroutine reset_CRMMScrb
