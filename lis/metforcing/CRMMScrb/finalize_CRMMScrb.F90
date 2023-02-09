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
! !MODULE: finalize_CRMMScrb
!  \label{finalize_CRMMScrb}
!
! !REVISION HISTORY:
!  5 Feb 2023: Sujay Kumar: Initial Implementation
! 
! !INTERFACE:
subroutine finalize_CRMMScrb(findex)
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use CRMMScrb_forcingMod, only : CRMMScrb_struc
!
! !DESCRIPTION:
!  Routine to cleanup CRMMS crb forcing related memory allocations.   
!
!  The arguments are: 
!  \begin{description}
!  \item[findex]
!    index of the forcing
!  \end{description}
!EOP
  implicit none
  
  integer   :: n
  integer   :: findex
  
  deallocate(CRMMScrb_struc)

end subroutine finalize_CRMMScrb
