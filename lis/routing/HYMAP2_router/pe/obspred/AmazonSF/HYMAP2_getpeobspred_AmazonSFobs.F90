!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: HYMAP2_getpeobspred_AmazonSFobs
!  \label{HYMAP2_getpeobspred_AmazonSFobs}
!
! !REVISION HISTORY:
! 02 Oct 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine HYMAP2_getpeobspred_AmazonSFobs(Obj_Func)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use HYMAP2_routingMod, only : HYMAP2_routing_struc
  use LIS_logMod,       only : LIS_verify, LIS_logunit

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: Obj_Func
!
! !DESCRIPTION:
!  
!  This routine retrieves the observation prediction, which is the
!  model's estimate of snow depth. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: sfField
  real, pointer          :: sf(:)
  integer                :: t,m
  integer                :: i
  integer                :: status


  n = 1

  call ESMF_StateGet(Obj_Func,"Amazon_SF",sfField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sfField,localDE=0,farrayPtr=sf,rc=status)
  call LIS_verify(status)

  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t = (i-1)*LIS_rc%nensem(n)+m 
        sf(t) = HYMAP2_routing_struc(n)%rivout(i,m)
     enddo
  enddo

end subroutine HYMAP2_getpeobspred_AmazonSFobs



