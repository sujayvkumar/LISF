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
! !ROUTINE: write_CNRS_em_obs_data
! \label{write_CNRS_em_obs_data}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine write_CNRS_em_obs_data(Obj_Space)
!!$! !USES: 
  use ESMF
!!$  use LIS_coreMod,  only : LIS_rc, LIS_masterproc
!!$  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
!!$       LIS_getNextUnitNumber, LIS_releaseUnitNumber
!!$  use LIS_fileIOMod,      only : LIS_create_output_directory
!!$  use LIS_historyMod,     only : LIS_writevar_gridded
!!$
!!$  implicit none
!!$! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!!$!
!!$! !DESCRIPTION:
!!$!  
!!$!  write the Walnut Gulch PBMR soil moisture data
!!$!  to disk
!!$!  
!!$!  The arguments are: 
!!$!  \begin{description}
!!$!  \item[n]    index of the nest
!!$!  \item[Obj\_State] observations state
!!$!  \end{description}
!!$!
!!$!EOP
!!$  real,    pointer    :: obsl(:)
!!$  type(ESMF_Field)    :: smField
!!$  type(ESMF_Logical)  :: data_update
!!$  character*100       :: obsname
!!$  integer             :: status
!!$  integer             :: ftn
!!$  integer             :: n 
!!$
!!$  n = 1
!!$
!!$  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
!!$       data_update, rc=status)
!!$  call LIS_verify(status)
!!$
!!$  if(data_update.eq..true.) then 
!!$     
!!$     call ESMF_StateGet(Obj_Space,"PBMR soil moisture",smField,&
!!$          rc=status)
!!$     call LIS_verify(status)
!!$
!!$     call ESMF_FieldGet(smField, localDE=0,farrayPtr=obsl,rc=status)
!!$     call LIS_verify(status)
!!$
!!$     if(LIS_masterproc) then 
!!$        ftn = LIS_getNextUnitNumber()
!!$        call CNRS_obsname(obsname)
!!$
!!$        call LIS_create_output_directory('PEOBS',style=2) 
!!$        open(ftn,file=trim(obsname), form='unformatted')
!!$     endif
!!$     
!!$     call LIS_writevar_gridded(ftn, n, obsl)
!!$     
!!$     if(LIS_masterproc) then 
!!$        call LIS_releaseUnitNumber(ftn)
!!$     endif
!!$
!!$  endif
!!$
end subroutine write_CNRS_em_obs_data

!BOP
! 
! !ROUTINE: CNRS_obsname
! \label{CNRS_obsname}
! 
! !INTERFACE: 
subroutine CNRS_obsname(obsname)
! !USES: 
  use LIS_coreMod, only : LIS_rc

! !ARGUMENTS: 
  character(len=*)     :: obsname
! 
! !DESCRIPTION: 
!  This method generates a timestamped filename for the processed
!  PBMR observations. 
! 
!EOP

  character(len=12)    :: cdate1

  write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
       LIS_rc%yr, LIS_rc%mo, &
       LIS_rc%da, LIS_rc%hr,LIS_rc%mn

  obsname = trim(LIS_rc%odir)//'/PEOBS/'//cdate1//'.1gs4r'
  
end subroutine CNRS_obsname



