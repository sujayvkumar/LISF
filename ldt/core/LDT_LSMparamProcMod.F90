!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
module LDT_LSMparamProcMod
!BOP
!
! !MODULE: LDT_LSMparamProcMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage the processing of 
!   LSM parameters
!
! !REVISION HISTORY: 
!  14 Aug 2014:  Sujay Kumar;  Initial Specification
! 
  use ESMF
  use LDT_coreMod
  use LDT_logMod

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_LSMparams_init
  public :: LDT_LSMparams_writeHeader
  public :: LDT_LSMparams_writeData

!BOP 
! 
! !ROUTINE: LDT_LSMparams_init
! \label{LDT_LSMparams_init}
! 
! !INTERFACE:
  interface LDT_LSMparams_init
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure LSMparams_init_LIS
     module procedure LSMparams_init_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing LSM parameters in both 
! in the preprocessing mode for LIS as well as in the LISHydro(WRFHydro) 
! preprocessing mode. 
!EOP 
  end interface

contains  
!BOP
! !ROUTINE: LSMparams_init_LIS
! \label{LSMparams_init_LIS}
!
! !INTERFACE: 
  subroutine LSMparams_init_LIS()

    integer :: flag
    integer :: k
    
    flag = 0

    do k=1,LDT_rc%nLSMs
       if(LDT_rc%lsm(k).ne."none") then          
          call lsmparamprocinit(trim(LDT_rc%lsm(k))//char(0),flag)
       endif
    enddo
    
  end subroutine LSMparams_init_LIS


!BOP
! !ROUTINE: LSMparams_init_LISHydro
! \label{LSMparams_init_LISHydro}
!
! !INTERFACE: 
  subroutine LSMparams_init_LISHydro(flag)
    
    integer   :: flag
    integer :: k

    flag = 1

    do k=1,LDT_rc%nLSMs
       if(LDT_rc%lsm(k).ne."none") then              
          call lsmparamprocinit(trim(LDT_rc%lsm(k))//char(0),flag)
       else
          write(LDT_logunit,*)'[ERR] Support for this LSM in the LISHydro preprocessing mode'
          write(LDT_logunit,*)'[ERR] is not supported'
          call LDT_endrun()
       endif
    enddo
    
  end subroutine LSMparams_init_LISHydro



!BOP
! !ROUTINE: LDT_LSMparams_writeHeader
! \label{LDT_LSMparams_writeHeader}
!
! !INTERFACE: 
  subroutine LDT_LSMparams_writeHeader(n,ftn,dimID, monthID)
    integer     :: n
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID
    integer     :: k
    
    do k=1,LDT_rc%nLSMs
       if(LDT_rc%lsm(k).ne."none") then              

          call lsmparamprocwriteheader(trim(LDT_rc%lsm(k))//char(0),&
               n,ftn,dimID, monthID)
       endif
    enddo
    
  end subroutine LDT_LSMparams_writeHeader


!BOP
! !ROUTINE: LDT_LSMparams_writeData
! \label{LDT_LSMparams_writeData}
!
! !INTERFACE: 
  subroutine LDT_LSMparams_writeData(n,ftn)

    integer     :: n
    integer     :: ftn
    integer     :: k

    do k=1,LDT_rc%nLSMs
       if(LDT_rc%lsm(k).ne."none") then              
          call lsmparamprocwritedata(trim(LDT_rc%lsm(k))//char(0),&
               n,ftn)
       endif
    enddo
    
  end subroutine LDT_LSMparams_writeData
     
end module LDT_LSMparamProcMod
