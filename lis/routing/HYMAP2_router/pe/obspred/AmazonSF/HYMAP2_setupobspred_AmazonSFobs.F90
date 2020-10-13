!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_setupobspred_AmazonSFobs
!  \label{HYMAP2_setupobspred_AmazonSFobs}
!
! !REVISION HISTORY:
! 2 Oct 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine HYMAP2_setupobspred_AmazonSFobs(OBSPred)
! !USES:
  use ESMF
  use HYMAP2_routingMod
  use LIS_coreMod     
  use LIS_logMod

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
!  This routine creates an entry in the Obs pred object from 
!  HYMAP2 used for parameter estimation
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: sfField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  sfField = ESMF_FieldCreate(arrayspec=realarrspec, &
       grid=LIS_vecRoutingTile(n), &
       name="Amazon_SF", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/sfField/),rc=status)
  call LIS_verify(status)

end subroutine HYMAP2_setupobspred_AmazonSFobs

