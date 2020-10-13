!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_plugins.h"
module LIS_routingoptue_pluginMod
!BOP
!
! !MODULE: LIS_routingoptue_pluginMod
!
! !DESCRIPTION:
!   This module contains the definition of the functions that
!   need to be defined to enable the use of a land surface
!   model in parameter estimation setup.
!
! !REVISION HISTORY:
!  16 Jul 09    Sujay Kumar  Initial Specification
!
  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LIS_routingoptue_plugin
!EOP

contains
!BOP
! !ROUTINE: LIS_routingoptue_plugin
!  \label{LIS_routingoptue_plugin}
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new ROUTING
! in a parameter estimation mode. The interface mandates that
! a number of routines be implemented and registered for
! each of the ROUTING that is used in a parameter estimation mode.
!
! !INTERFACE:
subroutine LIS_routingoptue_plugin
!EOP

#if ( ( defined OPTUE_ALG_ES )        || \
      ( defined OPTUE_ALG_LM )        || \
      ( defined OPTUE_ALG_GA )        || \
      ( defined OPTUE_ALG_SCEUA )     || \
      ( defined OPTUE_ALG_MCSIM )     || \
      ( defined OPTUE_ALG_RWMCMC )    || \
      ( defined OPTUE_ALG_DEMC )      || \
      ( defined OPTUE_ALG_DEMCZ ) )

   use LIS_pluginIndices
   use LIS_coreMod,  only : LIS_rc

#if ( defined ROUTE_HYMAP2_ROUTER)
   use HYMAP2_peMod, only : HYMAP2_setup_pedecvars
#endif

#if ( defined ROUTE_HYMAP2_ROUTER)

   external HYMAP2_set_pedecvars

   external HYMAP2_getpeobspred_AmazonSFobs
   external HYMAP2_setupobspred_AmazonSFobs

#endif


#if ( defined ROUTE_HYMAP2_ROUTER)

   call registerroutingpesetupdecisionspace(trim(LIS_hymap2routerid)//char(0), &
                                        HYMAP2_setup_pedecvars)
   call registerroutingpesetdecisionspace(trim(LIS_hymap2routerid)//char(0), &
                                      HYMAP2_set_pedecvars)

   call registerroutingpesetupobspred(trim(LIS_hymap2routerid)//"+"//      &
                                  trim(LIS_AmazonSFobsId)//char(0), &
                                  HYMAP2_setupobspred_AmazonSFobs)
   call registerroutingpegetobspred(trim(LIS_hymap2routerid)//"+"//      &
                                trim(LIS_AmazonSFobsId)//char(0), &
                                HYMAP2_getpeobspred_AmazonSFobs)
#endif
#endif
end subroutine LIS_routingoptue_plugin
end module LIS_routingoptue_pluginMod
