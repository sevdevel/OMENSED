SUBROUTINE Allocate_OMEN_Arrays
!************************************************************************
!
! *Allocate_OMEN_Arrays* Allocate memory for OMENSED sediment model arrays
!
! Author - Sebastiaan van de Velde
!
! Last update - 21 Jan 2021  @(OMENSED)Allocate_OMEN_Arrays.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMEN_model
!
! Module calls -
!
!************************************************************************

USE OMENpars
USE OMENvars

IMPLICIT NONE
!
!*Local variables
!
!
!1. OMENSED arrays
!------------------
!
ALLOCATE(OS_SWI_fluxes(1,MaxOMENSWIArids)) ! + M, + POC_tot, + POC_1, + POC_2
OS_SWI_fluxes = 0.0D0
ALLOCATE(OS_BW_conds(1,MaxOMENSWIArids))
OS_BW_conds = 0.0D0
!ALLOCATE(OS_part_fluxes(1,MaxOMENpartArids))
!OS_part_fluxes = 0.0D0
ALLOCATE(OS_boundaryconds(1,MaxOMENbcArids))
OS_boundaryconds = 0.0D0
ALLOCATE(OS_RCM_array(2,OS_RCM_fracs))
OS_RCM_array = 0.D0
!
!2. arrays for zTOC calculation
!------------------
!
ALLOCATE(OS_POC_conc_swi(OS_RCM_fracs))
OS_POC_conc_swi = 0.D0
ALLOCATE(aa1(OS_RCM_fracs))
aa1 = 0.D0
ALLOCATE(bb1(OS_RCM_fracs))
bb1 = 0.D0
ALLOCATE(aa2(OS_RCM_fracs))
aa2 = 0.D0
ALLOCATE(A11(OS_RCM_fracs))
A11 = 0.D0
ALLOCATE(A22(OS_RCM_fracs))
A22 = 0.D0
!
!3. arrays for diagenetic profiles
!------------------
!
!ALLOCATE(OS_POC_profile(OS_vertical_grid))
!OS_POC_profile = 0.D0
ALLOCATE(OS_z_vector(OS_vertical_grid))
OS_z_vector = 0.D0
ALLOCATE(OS_diag_profile(OS_vertical_grid,MaxOMENSWIArids)) 
OS_diag_profile=0.0D0

END SUBROUTINE Allocate_OMEN_Arrays

!========================================================================

SUBROUTINE Deallocate_OMEN_Arrays
!************************************************************************
!
! *Deallocate_OMEN_Arrays* Deallocate memory for OMENSED model arrays
!
! Author - Sebastiaan van de Velde
!
! Last update - 27 Jan 2021  @(OMENSED)Deallocate_OMEN_Arrays.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - CGEM_model
!
! Module calls -
!
!************************************************************************

USE OMENpars
USE OMENvars

IMPLICIT NONE

!
!1. OMENSED arrays
!------------------
!
DEALLOCATE(OS_SWI_fluxes) ! + M, + POC_tot, + POC_1, + POC_2
DEALLOCATE(OS_BW_conds)
!DEALLOCATE(OS_part_fluxes)
DEALLOCATE(OS_boundaryconds)
DEALLOCATE(OS_RCM_array)

!
!2. arrays for zTOC calculation
!------------------
!
DEALLOCATE(OS_POC_conc_swi)
DEALLOCATE(aa1)
DEALLOCATE(bb1)
DEALLOCATE(aa2)
DEALLOCATE(A11)
DEALLOCATE(A22)
!
!3. arrays for diagenetic profiles
!------------------
!
!DEALLOCATE(OS_POC_profile)
DEALLOCATE(OS_z_vector)
DEALLOCATE(OS_diag_profile) 


RETURN

END SUBROUTINE Deallocate_OMEN_Arrays
