!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah36_getsnodepvars
! \label{noah36_getsnodepvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!
! !INTERFACE:
subroutine noah36_getsnodepvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noah36_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField

  integer                :: t
  integer                :: status
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
 
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe(t) = noah36_struc(n)%noah(t)%sneqv
     snod(t) = noah36_struc(n)%noah(t)%snowh
  enddo
  
end subroutine noah36_getsnodepvars

