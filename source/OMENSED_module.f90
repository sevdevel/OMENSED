!************************************************************************
!  *OMENSED* Module containing the OMENSED diagenetic model (v1.0)
!
!  Author - Dominik Hülse and Sebastiaan van de Velde 
!
!  Version - @(OMENSED).f90  V0.1
!
! Description - UNITS: mol, cm, yr 
!
! Calling program -
! ******************************************************************************************************************************** !

MODULE OMENSED_module
    
USE OMENvars
USE OMENpars
USE OMENids
USE OMEN_Functions
USE OMEN_initialisation

IMPLICIT NONE

SAVE
    
CONTAINS

!========================================================================

SUBROUTINE OMENSED(SedVel, dum_sfcsumocn, &
                 & dum_new_swifluxes, dum_OMENSED_BC, dum_diag_profile, &
                 & dum_z_vector, dum_RCM_approx, dum_POC_conc_swi)

!************************************************************************
!
! *OMENSED* OMENSED model
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Last update - 28 Oct 2021  @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Main OMENSED subroutine: gets SWI POC wtpct and array of concentrations of solutes 
!               calls other subroutines and passes back fraction of POC preserved in sediments (dum_sed_pres_fracC)
!               and passess back calculated array of SWI fluxes of solutes
!
! Reference -
!
! Calling program - 
!
! Internal calls -
!
! Module calls -
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
!
REAL(8),INTENT(IN) :: SedVel                 ! sedimentation velocity (w) for OMEN-SED (cm yr-1)
REAL(8),DIMENSION(MaxOMENSWIArids),INTENT(IN) :: dum_sfcsumocn ! ocean composition interface array
REAL(8),DIMENSION(MaxOMENSWIArids),INTENT(INOUT) :: dum_new_swifluxes ! SWI return fluxes of solutes, 
                                                                   ! calculated with sediment-model 
                                                                   ! [pos values: flux from sediments to water-column]
REAL(8),DIMENSION(MaxOMENbcArids), INTENT(INOUT) :: dum_OMENSED_BC 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(INOUT) :: dum_RCM_approx
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(INOUT) :: dum_POC_conc_swi

!
!*Local variables
!
INTEGER :: i
REAL(8) :: loc_POC_flux_swi = 0.0D0    ! POC flux at SWI [mol/(cm^2 yr)]
REAL(8) :: loc_O2_swiflux    = 0.0D0    ! SWI return fluxes of O2 [mol/(cm^2 yr)]
REAL(8) :: loc_SO4_swiflux   = 0.0D0    ! SWI return fluxes of SO4 [mol/(cm^2 yr)]
REAL(8) :: loc_NO3_swiflux   = 0.0D0    ! SWI return fluxes of NO3 [mol/(cm^2 yr)]
REAL(8) :: loc_H2S_swiflux   = 0.0D0    ! SWI return fluxes of H2S [mol/(cm^2 yr)]
REAL(8) :: loc_NH4_swiflux   = 0.0D0    ! SWI return fluxes of H2S [mol/(cm^2 yr)]
REAL(8) :: loc_PO4_swiflux   = 0.0D0    ! SWI return fluxes of PO4 [mol/(cm^2 yr)]
REAL(8) :: loc_M_swiflux     = 0.0D0    ! SWI return fluxes of M - DOES NOT EXIST, JUST FOR DEBUGGING
REAL(8) :: loc_DIC_swiflux   = 0.0D0    ! SWI return fluxes of DIC [mol/(cm^2 yr)]
!REAL(8) :: loc_DIC_13C_swiflux = 0.0D0 ! SWI return fluxes of DIC [mol/(cm^2 yr)]
REAL(8) :: loc_ALK_swiflux   = 0.0D0    ! SWI return fluxes of ALK [mol/(cm^2 yr)]
REAL(8) :: loc_mixed_layer   = 5.0D0    ! mixed layer depth, to compare wt% with observations

REAL(8) :: dum_sed_mean_OM              ! mean OM wt% in upper mixed layer
REAL(8) :: dum_sed_OM_wtpc_bot          ! POC wt% at bottom of sediment column (geologic preservation)
REAL(8) :: dum_sed_pres_fracC           ! fraction POC-preserved/POC-deposited [-]
REAL(8) :: dum_sed_pres_fracP           ! fraction POP-preserved/POP-deposited [-]

! parameters for temperature dependent rate constants (as in John et al. 2014)
REAL(8) :: loc_T                        ! local temperature
REAL(8) :: loc_PO4_return_max
REAL(8) :: loc_fPOC                     ! Corg flux to the sediment [cm3 cm-2]
REAL(8) :: loc_fPOC_13C                 ! 13C of POC flux to the sediment [cm3 cm-2]

!==================================================================
! Give global arrays to local values
!==================================================================

loc_T = water_temp

dum_sed_OM_wtpc_bot = dum_OMENSED_BC(iarr_OS_POC_burial)
dum_sed_mean_OM  = dum_OMENSED_BC(iarr_OS_POC_swi)
dum_sed_pres_fracC = dum_OMENSED_BC(iarr_OS_POC_pres_frac)
dum_sed_pres_fracP = dum_OMENSED_BC(iarr_OS_POP_pres_frac)

dum_swiconc_O2  = dum_sfcsumocn(iarr_OS_O2)
dum_swiconc_NO3 = dum_sfcsumocn(iarr_OS_NO3)
dum_swiconc_NH4 = dum_sfcsumocn(iarr_OS_NH4)
dum_swiconc_SO4 = dum_sfcsumocn(iarr_OS_SO4)
dum_swiconc_H2S = dum_sfcsumocn(iarr_OS_H2S)
dum_swiconc_PO4 = dum_sfcsumocn(iarr_OS_PO4)
dum_swiconc_DIC = dum_sfcsumocn(iarr_OS_DIC)
dum_swiconc_ALK = dum_sfcsumocn(iarr_OS_AT)

loc_fPOC = dum_sfcsumocn(iarr_OS_POC)
loc_POC_flux_swi = conv_POC_cm3_mol*loc_fPOC
!print*, 'loc_FPOC= ',loc_POC_flux_swi
w = SedVel  
   
!==================================================================
! Check for low sedimentation velocities
!==================================================================

! Model crashed for low sediment accumulation rates, therefore:
! low sedimentation velocity -> Remineralize everything manually

IF((w .LE. const_real_nullsmall) .OR. OS_remineralize_all) THEN
    !!! Remineralize everything manually
    dum_sed_pres_fracC = 0.0D0        ! sed POC preservation to zero
    dum_sed_pres_fracP = 0.0D0
    dum_sed_OM_wtpc_bot = 0.0D0
    loc_O2_swiflux = conv_POC_cm3_mol*loc_fPOC*(-OC/SD)
    loc_NO3_swiflux = 0.0D0
    loc_NH4_swiflux = 0.0D0
    loc_SO4_swiflux = 0.0D0 !conv_POC_cm3_mol*loc_fPOC*(-SO4C/SD)
    loc_H2S_swiflux = -loc_SO4_swiflux
    loc_PO4_swiflux = conv_POC_cm3_mol*loc_fPOC*PC1/SD
    loc_DIC_swiflux = conv_POC_cm3_mol*loc_fPOC*DICC1/SD
    loc_ALK_swiflux = 2.0D0*loc_H2S_swiflux + conv_POC_cm3_mol*loc_fPOC*ALKROX/SD   !16/106

    GOTO 100

END IF 

! Model crashed for low sediment accumulation rates, therefore:
! too low sedimentation velocity -> set to minimum

IF (w .LE. 5.0D-4) THEN
    w = 5.0D-4     !(5.0e-4 for OAE2; 4.0e-4 for modern)
END IF

CALL initialise_OMEN_params(WaterTemp=loc_T, O2_BW=dum_swiconc_O2, sed_w=w)

IF (nt.EQ.0) GOTO 200

!==================================================================
! Deconvolute RCM into multifractions
!==================================================================

CALL RCM_approx(RCMarray=dum_RCM_approx,RCM_a=RCM_a,RCM_nu=RCM_nu)

!DO i=1, OS_RCM_fracs
!   print*, 'F', OS_RCM_array(1,i)
!   print*, 'k', OS_RCM_array(2,i)
!END DO
!print*, 'Fsum', SUM(OS_RCM_array(1,:))

 ! Check for no POC deposited -> nothing preserved
IF (loc_POC_flux_swi .LE. const_real_nullsmall) THEN
    
    dum_sed_pres_fracC = 0.0D0
    dum_sed_pres_fracP = 0.0D0
    dum_sed_OM_wtpc_bot = 0.0D0
    loc_O2_swiflux = 0.0D0
    loc_NO3_swiflux = 0.0D0
    loc_SO4_swiflux = 0.0D0
    loc_NH4_swiflux = 0.0D0
    loc_H2S_swiflux = 0.0D0
    loc_PO4_swiflux = 0.0D0
    loc_DIC_swiflux = 0.0D0
    loc_ALK_swiflux = 0.0D0
         
    GOTO 100
END IF
          
CALL OMENSED_calculate_zTOC(dum_POC_flux_swi=loc_POC_flux_swi,dum_sed_pres_fracC=dum_sed_pres_fracC, &
                          & dum_sed_OM_wtpc_bot=dum_sed_OM_wtpc_bot,dum_POC_conc_swi=dum_POC_conc_swi, &
                          & dum_diag_profile=dum_diag_profile, dum_z_vector=dum_z_vector, &
						  & dum_RCM_approx=dum_RCM_approx)

! CHECK IF TOC preservation results in insane values, i.e. everything remineralized
! Then calculate SWI-fluxes "manually"
IF ((dum_sed_pres_fracC.NE.dum_sed_pres_fracC).OR.(dum_sed_pres_fracC.LE.0.0D0) &
   & .OR.(dum_sed_pres_fracC .GT. 1.0D0)) THEN
    print*, 'TOO FAST', dum_sed_pres_fracC   
        
    dum_sed_pres_fracC = 0.0D0        ! sed TOC preservation to zero
    dum_sed_pres_fracP = 0.0D0
    dum_sed_OM_wtpc_bot = 0.0D0
    loc_O2_swiflux = conv_POC_cm3_mol*loc_fPOC*(-OC/SD)
    loc_NO3_swiflux = 0.0D0
    loc_NH4_swiflux = 0.0D0
    loc_SO4_swiflux = 0.0D0 !conv_POC_cm3_mol*loc_fPOC*(-SO4C/SD)
    loc_H2S_swiflux = -loc_SO4_swiflux
    loc_PO4_swiflux = conv_POC_cm3_mol*loc_fPOC*PC1/SD
    loc_DIC_swiflux = conv_POC_cm3_mol*loc_fPOC*DICC1/SD
    loc_ALK_swiflux = 2.0D0*loc_H2S_swiflux + conv_POC_cm3_mol*loc_fPOC*ALKROX/SD   !16/106
    
    GOTO 100
    
END IF 

IF (dum_swiconc_O2 .LE. const_real_nullsmall) THEN
    loc_O2_swiflux = 0.0D0            ! if negative [O2] -> no SWI flux
ELSE
    print*, 'BW_O2 ztoc', OS_BW_conds(1,iarr_OS_O2)   
    CALL OMENSED_calculate_zO2(dum_swiconc_O2=dum_swiconc_O2, dum_swiflux_O2=loc_O2_swiflux, &
                             & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                             & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
    IF (loc_O2_swiflux .GE. 0.0D0) THEN
        print*,'---------- loc_O2_swiflux positiv ----------', loc_O2_swiflux
        loc_O2_swiflux = 0.0D0
    END IF
END IF

IF (dum_swiconc_NO3.GT.const_real_nullsmall) THEN
    CALL OMENSED_calculate_zNO3(dum_swiconc_NO3=dum_swiconc_NO3, dum_swiflux_NO3=loc_NO3_swiflux, &
                              & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                              & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
ELSE
    loc_NO3_swiflux = 0.0D0  
END IF
    
IF (dum_swiconc_SO4.GT.const_real_nullsmall) THEN
    CALL OMENSED_calculate_zSO4(dum_swiconc_SO4=dum_swiconc_SO4, dum_swiflux_SO4=loc_SO4_swiflux, &
                              & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                              & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
    IF (loc_SO4_swiflux.GT.0.0D0) THEN
        print*,'---------- loc_SO4_swiflux positiv ----------', loc_SO4_swiflux
        loc_SO4_swiflux = 0.0D0
    END IF
ELSE
    zso4 = zno3
END IF
             
CALL OMENSED_calculate_zNH4(dum_swiconc_NH4=dum_swiconc_NH4, dum_swiflux_NH4=loc_NH4_swiflux, &
                          & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                              & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
             
CALL OMENSED_calculate_zH2S(dum_swiconc_H2S=dum_swiconc_H2S, dum_swiflux_H2S=loc_H2S_swiflux, &
                          & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                              & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
              
! calculate maximum return flux of PO4 from Porg in sediment cell
loc_PO4_return_max = loc_fPOC*conv_POC_cm3_mol*PC1/SD
             
IF (par_sed_huelse2017_P_cycle) THEN
    ! explicit PO4 calculation
    CALL OMENSED_calculate_zPO4_M(dum_swiconc_PO4=dum_swiconc_PO4, dum_swiflux_PO4=loc_PO4_swiflux, &
                                & dum_swiflux_M=loc_M_swiflux, &
                                & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                                & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
    dum_sed_pres_fracP = (loc_PO4_return_max-loc_PO4_swiflux)/loc_PO4_return_max
ELSE
    IF (par_sed_huelse2017_sim_P_loss) THEN
        ! PO4 hack: remineralise just the non-sulfurised POC and calculate PO4 return flux
        loc_PO4_swiflux = loc_POC_flux_swi*PC1/SD
        dum_sed_pres_fracP = (loc_PO4_return_max-loc_PO4_swiflux)/loc_PO4_return_max
        IF (par_sed_huelse2017_sim_P_loss_pres_fracC) THEN        ! preserve same fraction as Corg
            loc_PO4_swiflux = loc_fPOC*conv_POC_cm3_mol*PC1/SD*(1-dum_sed_pres_fracC)
            dum_sed_pres_fracP = dum_sed_pres_fracC
        END IF
    ! simulate P-regeneration under anoxia:
    ELSE                                                                
        ! PO4 hack: remineralise all POC and calculate PO4 return flux
        loc_PO4_swiflux = loc_fPOC*conv_POC_cm3_mol*PC1/SD
        dum_sed_pres_fracP = 0.0D0
    END IF
END IF  ! par_sed_huelse2017_P_cycle                       

CALL OMENSED_calculate_zDIC(dum_swiconc_DIC=dum_swiconc_DIC, dum_swiflux_DIC=loc_DIC_swiflux, &
                          & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                          & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
                
CALL OMENSED_calculate_zALK(dum_swiconc_ALK=dum_swiconc_ALK, dum_swiflux_ALK=loc_ALK_swiflux, &
                          & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx, &
                          & dum_z_vector=dum_z_vector, dum_diag_profile=dum_diag_profile)
                              
loc_ALK_swiflux = loc_ALK_swiflux + dum_POC_total_flux_zinf*ALKROX/SD

100 CONTINUE

! calculate mean OM concentration [wt%] in upper x cm
! convert from mol/cm3 to wt%
! just two fractions:
dum_sed_mean_OM = 1.0D0/loc_mixed_layer * 100.0D0*12.0D0/rho_sed* &
     & FUN_OMENSED_calcOM(zU=0.0D0, zL=loc_mixed_layer, reac=1.0D0, dum_POC_conc_swi=dum_POC_conc_swi)

!==================================================================
! Give local values back to global arrays
!==================================================================

! Fluxes
dum_new_swifluxes(iarr_OS_POC)  = loc_POC_flux_swi!*1.D3*1.D4/365.25D0      
  
dum_new_swifluxes(iarr_OS_O2)  = loc_O2_swiflux!*1.D3*1.D4/365.25D0        
dum_new_swifluxes(iarr_OS_NO3) = loc_NO3_swiflux!*1.D3*1.D4/365.25D0     
dum_new_swifluxes(iarr_OS_NH4) = loc_NH4_swiflux!*1.D3*1.D4/365.25D0  
dum_new_swifluxes(iarr_OS_SO4) = loc_SO4_swiflux!*1.D3*1.D4/365.25D0    
dum_new_swifluxes(iarr_OS_H2S) = loc_H2S_swiflux!*1.D3*1.D4/365.25D0   
dum_new_swifluxes(iarr_OS_PO4) = loc_PO4_swiflux!*1.D3*1.D4/365.25D0 
dum_new_swifluxes(iarr_OS_DIC) = loc_DIC_swiflux!*1.D3*1.D4/365.25D0  
dum_new_swifluxes(iarr_OS_AT)  = loc_ALK_swiflux!*1.D3*1.D4/365.25D0   

! Boundary conditions
dum_OMENSED_BC(iarr_OS_w)             = w
dum_OMENSED_BC(iarr_OS_zox)           = zox
dum_OMENSED_BC(iarr_OS_zno3)          = zno3
dum_OMENSED_BC(iarr_OS_zSO4)          = zSO4
dum_OMENSED_BC(iarr_OS_POC_burial)    = dum_sed_OM_wtpc_bot
dum_OMENSED_BC(iarr_OS_POC_swi)       = dum_sed_mean_OM
dum_OMENSED_BC(iarr_OS_POC_pres_frac) = dum_sed_pres_fracC
dum_OMENSED_BC(iarr_OS_POP_pres_frac) = dum_sed_pres_fracP

IF (print_results) THEN
   print*,'Temp C =', dum_sfcsumocn(iarr_OS_T) - 273.15D0
   print*,'loc_POC_flux_swi = ', loc_POC_flux_swi*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'Fraction POC-pres =' , dum_sed_pres_fracC
   print*,'(1-presC)*POCin   = ', (1-dum_sed_pres_fracC)*(loc_POC_flux_swi) &
                                & *1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL DIC SWI flx = ', dum_new_swifluxes(iarr_OS_DIC)*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'loc_sed_burial    = ', w
   print*,'dum_swiconc_O2 = ', dum_swiconc_O2
   print*,'dum_swiconc_SO4 = ', dum_swiconc_SO4
   print*,'dum_swiconc_H2S = ', dum_swiconc_H2S
   print*,'dum_swiconc_PO4 = ', dum_swiconc_PO4
   print*,'dum_swiflux_M = ', dum_swiflux_M
   print*,'dum_swiconc_DIC = ', dum_swiconc_DIC
   print*,'dum_swiconc_ALK = ', dum_swiconc_ALK
   print*,' '
   print*,'zox = ', zox
   print*,'FINAL O2 SWI flux = ', dum_new_swifluxes(iarr_OS_O2)*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'zno3 = ', zno3
   print*,'FINAL NO3 SWI flux = ', loc_NO3_swiflux*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'zso4 = ', zso4
   print*,'FINAL SO4 SWI flux = ', dum_new_swifluxes(iarr_OS_SO4)*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL NH4 SWI flux = ', loc_NH4_swiflux*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL H2S SWI flux = ', dum_new_swifluxes(iarr_OS_H2S)*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL PO4 SWI flux = ', dum_new_swifluxes(iarr_OS_PO4)*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL M SWI flux = ', loc_M_swiflux*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL DIC SWI flux = ', dum_new_swifluxes(iarr_OS_DIC)*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'FINAL ALK SWI flux = ', loc_ALK_swiflux*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
   print*,'Fraction POC-preserved/POC-deposited =' , dum_sed_pres_fracC
   print*,' '

END IF

CALL Conservation_Check(fluxarray=dum_new_swifluxes,bcarray=dum_OMENSED_BC)


200 CONTINUE

RETURN

END SUBROUTINE OMENSED
    
!========================================================================

SUBROUTINE OMENSED_calculate_zTOC(dum_POC_flux_swi, dum_sed_pres_fracC, & ! MULTIG
                                & dum_sed_OM_wtpc_bot,dum_POC_conc_swi, &
                                & dum_diag_profile, dum_z_vector, dum_RCM_approx) 
!************************************************************************
!
! *OMENSED_calculate_zTOC* calculate benthic burial/recycling fluxes 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
!       
REAL(8), INTENT(IN)    :: dum_POC_flux_swi    ! POC flux at SWI   [mol/(cm2 yr)]
REAL(8), INTENT(INOUT) :: dum_sed_pres_fracC  ! POC concentrations at zinf
REAL(8), INTENT(INOUT) :: dum_sed_OM_wtpc_bot ! POC wt% at bottom of sediment column (geologic preservation)
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(INOUT) :: dum_POC_conc_swi ! concentration of POC fractions at SWI
!REAL(8), DIMENSION(OS_vertical_grid), INTENT(INOUT) :: dum_POC_profile
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8), DIMENSION(OS_RCM_fracs) :: loc_POC_conc_zinf
REAL(8), DIMENSION(OS_RCM_fracs) :: dC1dz, C1flx,Cflx, F_TOC1         ! Cflx: Sed input flux to upper boundary, per cm^2 water column
REAL(8) :: F_TOC, loc_POC_C0                                                 ! Flux through lower boundary zinf, per cm^2 water-column
REAL(8), DIMENSION(OS_RCM_fracs) :: loc_POC_conc_swi_nonbio, loc_POC_conc_swi, loc_POC_flux
REAL(8), DIMENSION(OS_RCM_fracs,OS_vertical_grid) :: loc_POC_profile
INTEGER :: i, j

print*,' ------------------ START zTOC ---------------------'

! calculate concentration (mol/cm3) at very top of sediments not accounting for biodiffusion (just account for advection)
! NOTE: this is needed when calculating final fraction of POC preservation
! if taking [POC] at SWI including biodiffusion term (e.g. dum_POC1_conc_swi) 
! the preservation is too high (as not all incoming POC is accounted for)

loc_POC_C0 = dum_POC_flux_swi*1.0D0/(1.0D0-por)*1.0D0/w
!print*, 'loc_POC_C0 =', loc_POC_C0

DO i=1, OS_RCM_fracs
   
   loc_POC_conc_swi_nonbio(i) = dum_RCM_approx(1,i)*loc_POC_C0
   !print*, 'loc_POC_conc_swi_nonbio', i, ' =', loc_POC_conc_swi_nonbio(i)
   
   aa1(i) = (w-DSQRT(w**2.0D0+4.0D0*DC1*dum_RCM_approx(2,i)))/(2.0D0*DC1)
   bb1(i) = (w+DSQRT(w**2.0D0+4.0D0*DC1*dum_RCM_approx(2,i)))/(2.0D0*DC1)
   aa2(i) = (-OS_RCM_array(2,i)/w)
   !print*, 'aa1', i, ' =', aa1(i)
   !print*, 'bb1', i, ' =', bb1(i)
   !print*, 'aa2', i, ' =', aa2(i)
   
   loc_POC_flux(i)=dum_RCM_approx(1,i)*dum_POC_flux_swi
   loc_POC_conc_swi(i) = (loc_POC_flux(i)*(-aa1(i)*DEXP(aa1(i)*zbio)+bb1(i)*DEXP(bb1(i)*zbio)))/ &
                  & (-DC1*bb1(i)*aa1(i)*DEXP(bb1(i)*zbio) + DC1*bb1(i)*aa1(i)*DEXP(aa1(i)*zbio) + &
                  & DC1*bb1(i)*aa1(i)*por*DEXP(bb1(i)*zbio) - DC1*bb1(i)*aa1(i)*por*DEXP(aa1(i)*zbio) - &
                  & w*aa1(i)*DEXP(aa1(i)*zbio) + &
                  & w*bb1(i)*DEXP(bb1(i)*zbio) + w*por*aa1(i)*DEXP(aa1(i)*zbio) - w*por*bb1(i)*DEXP(bb1(i)*zbio))
   
   !print*, 'loc_POC_conc_swi', i, ' =', loc_POC_conc_swi(i)
   ! calculate TOC SWI concentration from flux
   A11(i) = -(loc_POC_conc_swi(i)*bb1(i)*DEXP(bb1(i)*zbio))/(aa1(i)*DEXP(aa1(i)*zbio)-bb1(i)*DEXP(bb1(i)*zbio)+OS_tolcte) 
   A22(i)=(A11(i)*(DEXP(aa1(i)*zbio)-DEXP(bb1(i)*zbio))+loc_POC_conc_swi(i)*DEXP(bb1(i)*zbio))/(DEXP(aa2(i)*zbio)+OS_tolcte) 
   !print*, 'A11', i, ' =', A11(i)
   !print*, 'A22', i, ' =', A22(i)
   
   ! Cflx: Sed input flux to upper boundary, per cm^2 water column
   IF (z0.LT.zbio) THEN
      dC1dz(i) =  A11(i)*(aa1(i)*DEXP(aa1(i)*z0)-bb1(i)*DEXP(bb1(i)*z0))+loc_POC_conc_swi(i)*bb1(i)*DEXP(bb1(i)*z0)
      C1flx(i) = - (1.0D0-por)*(-DC1*dC1dz(i) + w*loc_POC_conc_swi(i))
   ELSE
      C1flx(i) = - (1.0D0-por)*w*loc_POC_conc_swi(i)
   END IF 
   
   !print*, 'C1flx', i, ' =', C1flx(i)
   Cflx = SUM(C1flx)
   
   ! Flux through lower boundary zinf, per cm^2 water-column 
   F_TOC1(i) = (1-por)*w*A22(i)*DEXP(aa2(i)*zinf)
   F_TOC = SUM(F_TOC1)
   
   ! Concentration at lower boundary zinf
   IF (zinf.LT.zbio) THEN
      loc_POC_conc_zinf(i)=A11(i)*(DEXP(aa1(i)*zinf)-DEXP(bb1(i)*zinf))+loc_POC_conc_swi(i)*DEXP(bb1(i)*zinf)
   ELSE
      loc_POC_conc_zinf(i)=A22(i)*DEXP(aa2(i)*zinf)
   END IF
   
   DO j=2, OS_vertical_grid-1
      IF (dum_z_vector(j).LE.zbio) THEN
         loc_POC_profile(i,j) = A11(i)*(DEXP(aa1(i)*dum_z_vector(j))-DEXP(bb1(i)*dum_z_vector(j)))+ &
                              & loc_POC_conc_swi(i)*DEXP(bb1(i)*dum_z_vector(j))
      ELSE
         loc_POC_profile(i,j) = A22(i)*DEXP(aa2(i)*dum_z_vector(j))
      END IF
   END DO

   dum_POC_conc_swi(i) = loc_POC_conc_swi(i)

END DO

! DH: need to give back fraction buried of initially deposited (so fraction of the input values to this subroutine)
dum_diag_profile(1,iarr_OS_POC)                = SUM(loc_POC_conc_swi)!SUM(loc_POC_conc_swi_nonbio)
dum_diag_profile(OS_vertical_grid,iarr_OS_POC) = SUM(loc_POC_conc_zinf)
DO j=2, OS_vertical_grid-1
   dum_diag_profile(j,iarr_OS_POC) = SUM(loc_POC_profile(:,j))
   !print*, 'zTOC= ',dum_diag_profile(j,iarr_OS_POC)
   !print*, 'j= ',j
END DO

dum_sed_pres_fracC = SUM(loc_POC_conc_zinf)/SUM(loc_POC_conc_swi_nonbio) 
dum_POC_total_flux_zinf = F_TOC
dum_sed_OM_wtpc_bot = 100.0D0*12.0D0/rho_sed*SUM(loc_POC_conc_zinf)

RETURN

END SUBROUTINE OMENSED_calculate_zTOC

!========================================================================

SUBROUTINE OMENSED_calculate_zO2(dum_swiconc_O2, dum_swiflux_O2, dum_POC_conc_swi, dum_RCM_approx, &
                               & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zO2* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
!       
REAL(8),INTENT(IN) :: dum_swiconc_O2          ! O2 concentrations at SWI
REAL(8),INTENT(INOUT) :: dum_swiflux_O2   ! O2 flux: TODO check! (+) sediment -> bottom waters
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
!
REAL(8) :: flxzox, conczox, loc_conczinf, fun0, zL, tol
INTEGER :: bctype

print*, '---------------------- START zO2 ------------------------ '

zox = 1D-10
flxzox = 0.0D0
conczox = 0.0D0
loc_conczinf = 0.0D0
dum_swiflux_O2 = 0.0D0
fun0 = 0.0D0
zL = 0.0D0
tol = 0.0D0

fun0 = FUN_OMENSED_zO2(z=zox)

! Try zero flux at zinf and see if we have any O2 left
bctype = 2

CALL OMENSED_calculate_zO2_calcbc(zox=zinf, bctype=bctype, flxzox=flxzox, conczox=loc_conczinf, &
                                & flxswi=dum_swiflux_O2, r_zxf=r_zxf, dum_POC_conc_swi=dum_POC_conc_swi, &
                                & dum_RCM_approx=dum_RCM_approx, dum_z_vector=dum_z_vector, &
                                & dum_diag_profile=dum_diag_profile)

IF (fun0 .GE. 0.0D0) THEN   ! eg zero oxygen at swi
   zox = 0.0D0   ! DH 241016 was 1e-10
   bctype = 1
   loc_conczinf = 0.0D0
ELSEIF (loc_conczinf .GE. 0.0D0) THEN      ! still O2 at zinf -> zox = zinf
   zox = zinf
   zno3 = zinf
   bctype = 2
ELSE                        ! search zox in the interval
   bctype = 1
   zL=1D-10
   tol=1D-16
   zox = FUN_zbrent(FUN_OMENSED_zO2, zL, zinf, tol) ! SVDV: FUN_zbrent(fun0, zL, zinf, tol) ???
   zno3 = zox
   loc_conczinf = 0.0D0
END IF

CALL OMENSED_calculate_zO2_calcbc(zox=zox, bctype=bctype, flxzox=flxzox, conczox=conczox, &
                                & flxswi=dum_swiflux_O2, r_zxf=r_zxf, dum_POC_conc_swi=dum_POC_conc_swi, &
                                & dum_RCM_approx=dum_RCM_approx, dum_z_vector=dum_z_vector, &
                                & dum_diag_profile=dum_diag_profile)
                                      
dum_swiflux_O2 = dum_swiflux_O2 - por*w*(dum_swiconc_O2-loc_conczinf)

END SUBROUTINE OMENSED_calculate_zO2

!========================================================================

SUBROUTINE OMENSED_calculate_zO2_calcbc(zox, bctype, flxzox, conczox, flxswi,r_zxf, &
                                      & dum_POC_conc_swi, dum_RCM_approx, &
                                      & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zO2_calcbc* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED_calculate_zO2 / FUN_zO2
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
! 
REAL(8), INTENT(IN) :: zox
INTEGER, INTENT(IN) :: bctype
REAL(8), INTENT(INOUT) :: flxzox, conczox, flxswi, r_zxf
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
!
REAL(8) :: reac !ls , z0, zox
INTEGER :: ltype
REAL(8) :: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
REAL(8) :: e_0, dedz_0, f_0, dfdz_0, g_0, dgdz_0
REAL(8) :: e_zox, dedz_zox, f_zox, dfdz_zox, g_zox, dgdz_zox
REAL(8) :: rO2_AO2, rO2_BO2, Dzox, Dswi
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
INTEGER :: j

!   reactive terms: OM degradation (-) and nitrification (-)
reac = -OC-2.0D0*gammaNH4*NC1

! calculate solution for given zox
! Preparation: sort out solution-matching across bioturbation boundary (if necessary)
CALL sub_prepfg_l12(reac=reac, ktemp=0.0D0, zU=z0, zL=zox, D1=DO21, D2=DO22, ls_a=ls_a, ls_b=ls_b, ls_c=ls_c, &
                  & ls_d=ls_d, ls_e=ls_e, ls_f=ls_f, ltype=ltype, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! basis functions at upper boundary
!        print*, 'z0, reac1, reac2 ', z0, reac1, reac2
!        print*, 'ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, DO21, DO22', ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, DO21, DO22
!        print*, 'ltype, DO21, DO22, ltype', ltype, DO21, DO22, ltype
CALL sub_calcfg_l12(z=z0, reac=reac, ktemp=0.0D0, ls_a=ls_a, ls_b=ls_b, ls_c=ls_c, ls_d=ls_d, ls_e=ls_e, ls_f=ls_f, &
                  & ls_D1=DO21, ls_D2=DO22, ltype=ltype, e=e_0, dedz=dedz_0, f=f_0, dfdz=dfdz_0, &
                  & g=g_0, dgdz=dgdz_0, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                  
! ... and lower boundary
CALL sub_calcfg_l12(z=zox, reac=reac, ktemp=0.0D0, ls_a=ls_a, ls_b=ls_b, ls_c=ls_c, ls_d=ls_d, ls_e=ls_e, ls_f=ls_f, &
                  & ls_D1=DO21, ls_D2=DO22, ltype=ltype, e=e_zox, dedz=dedz_zox, f=f_zox, dfdz=dfdz_zox, &
                  & g=g_zox, dgdz=dgdz_zox, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)


! Solve for AO2, BO2 given boundary conditions (expressed in terms of transformed soln)

IF (bctype.EQ.1) THEN
   ! Case 1 zero concentration at zox
   ! AO2*e_zox   +   BO2*f_zox  + g_zox = 0;
   ! AO2*e_0     +   BO2*f_0     + g_0  = swi.O20;

   ! | e_zox f_zox |  |AO2|   = | -g_zox         |
   ! | e_0   f_0   |  |BO2|     | swi.O20 - gt_0 |

   CALL sub_solve2eqn(e_zox, f_zox, e_0, f_0, -g_zox, dum_swiconc_O2 - g_0, rO2_AO2, rO2_BO2)

ELSE
   ! Case  2 zero flux at zox
   ! AO2*dedz_zox +  BO2*dfz_zox + dgz_zox = 0;
   ! AO2*e_0     +   BO2*f_0     + g_0     = swi.O20;
   ! a            b       c   d       e           f
   !            print*, 'dedz_zox, dfdz_zox, e_0, f_0, -dgdz_zox, dum_swiconc_O2 - g_0:', dedz_zox, dfdz_zox, e_0, f_0, -dgdz_zox, dum_swiconc_O2 - g_0
   !            print*, 'dum_swiconc_O2, g_0:', dum_swiconc_O2, g_0
   CALL sub_solve2eqn(dedz_zox, dfdz_zox, e_0, f_0, -dgdz_zox, dum_swiconc_O2 - g_0, rO2_AO2, rO2_BO2)

END IF

IF (zox .LT. zbio) THEN
   Dzox = DO21
ELSE
   Dzox = DO22
END IF

flxzox =  Dzox*(rO2_AO2*dedz_zox + rO2_BO2*dfdz_zox + dgdz_zox)         ! no por factor as this is per cm^2 pore area
conczox = rO2_AO2*e_zox + rO2_BO2 * f_zox + g_zox

IF (0.D0.LT.zbio) THEN
   Dswi = DO21
ELSE
   Dswi = DO22
END IF

flxswi = por*(Dswi*(rO2_AO2*dedz_0+rO2_BO2*dfdz_0 + dgdz_0)) ! just diffusive flux - w*dum_swiconc_O2)   ! por fac so this is per cm^2 water column
r_zxf = zox/(zoxgf + zox) ! + const_real_nullsmall)   ! roll off oxidation at low zox


! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_O2) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zox) THEN
   
      CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac, ktemp=0.0D0, ls_a=ls_a, ls_b=ls_b, ls_c=ls_c, ls_d=ls_d, &
                        & ls_e=ls_e, ls_f=ls_f, ls_D1=DO21, ls_D2=DO22, ltype=ltype, &
                        & e=e_zx, dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                        & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
      dum_diag_profile(j,iarr_OS_O2) = rO2_AO2*e_zx + rO2_BO2 * f_zx + g_zx
    END IF
    
END DO


IF ((flxzox.NE.flxzox).OR.(conczox.NE.conczox).OR.(flxswi.NE.flxswi)) THEN !check for NaN if then give value as in matlab.....
    print*,' '
    print*,' '
    print*,'------ zO2_calcbc --------- flxzox is INFFFFFFFFFFFFFFFF'
    print*,' zox tested', zox
    print*,'flxzox ', flxzox
    print*,'conczox ', conczox
    print*,'flxswi ', flxswi
    print*,'dum_swiconc_O2 ', dum_swiconc_O2
    print*,' '

    conczox = -1.0D0
    flxzox = -1.0D0
END IF

END SUBROUTINE OMENSED_calculate_zO2_calcbc

!========================================================================

FUNCTION FUN_OMENSED_calcFO2(z)
!************************************************************************
!
! *FUN_OMENSED_calcFO2* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_zO2
!
!************************************************************************
         
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: z
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcFO2, tmpreac

FUN_OMENSED_calcFO2 = 0.0D0
! Oxidation of reduced species at zox (NEED A RATIO for ODU! and add NH4 adsorption!)
tmpreac=gammaH2S*O2H2S*SO4C+2.0D0*gammaNH4*NC1

FUN_OMENSED_calcFO2 = z/(zoxgf + z) * FUN_OMENSED_calcReac(zU=z, zL=zinf, reac=tmpreac, dum_RCM_approx=OS_RCM_array, &
                                                           & dum_POC_conc_swi=OS_POC_conc_swi)   ! had in denominator: ... + const_real_nullsmall

END FUNCTION FUN_OMENSED_calcFO2

!========================================================================

FUNCTION FUN_OMENSED_zO2(z)
!************************************************************************
!
! *FUN_OMENSED_calcFO2* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_zO2
!
!************************************************************************
         
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: z
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_zO2, flxzox, conczox, flxswi, r_zxf
INTEGER :: bctype
        
flxzox = 0.0D0
conczox = 0.0D0
flxswi = 0.0D0

bctype = 1

CALL OMENSED_calculate_zO2_calcbc(zox=z, bctype=bctype, flxzox=flxzox, conczox=conczox, &
                                & flxswi=flxswi, r_zxf=r_zxf, dum_POC_conc_swi=OS_POC_conc_swi, &
                                & dum_RCM_approx=OS_RCM_array, dum_z_vector=OS_z_vector, &
                                & dum_diag_profile=OS_diag_profile)

FUN_OMENSED_zO2 = flxzox + FUN_OMENSED_calcFO2(z)

END FUNCTION FUN_OMENSED_zO2

!========================================================================

SUBROUTINE OMENSED_calculate_zNO3(dum_swiconc_NO3, dum_swiflux_NO3, dum_POC_conc_swi, dum_RCM_approx, &
                                & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zNO3* calculate nitrate penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
! 
REAL(8), INTENT(IN) :: dum_swiconc_NO3                ! NO3 concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_NO3         ! NO3 flux: TODO check! (+) sediment -> bottom waters
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAL(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
!
REAL(8) :: flxzinf, conczinf, flxzno3, conczno3, zL, tol
INTEGER :: bctype

print*, '---------------------- START zNO3 ------------------------------- '

! Try zero flux at zinf and see if we have any NO3 left - in
! the rare case of zox < zinf but zNO3 = zinf, also
! calculate [NO3] at zinf for advective loss
bctype = 2

CALL OMENSED_calculate_zNO3_calcbc(zNO3=zinf, bctype=bctype, flxzno3=flxzinf, conczno3=conczinf, &
                                 & dum_swiflux_NO3=dum_swiflux_NO3, dum_POC_conc_swi=dum_POC_conc_swi, &
                                 & dum_RCM_approx=dum_RCM_approx, dum_z_vector=dum_z_vector, &
                                 & dum_diag_profile=dum_diag_profile)

IF (zox.EQ.zinf) THEN
   zno3 = zinf
   bctype = 2
ELSE
   IF (conczinf.GT.0.0D0) THEN
      zno3 = zinf
      bctype = 2;
   ELSE
      zL=1D-10
      tol=1D-16
      zno3 = FUN_zbrent(FUN_OMENSED_zNO3, max(zL,zox), zinf, tol)
      bctype = 1;
      conczinf = 0.0D0
   END IF
END IF 

CALL OMENSED_calculate_zNO3_calcbc(zNO3=zno3, bctype=bctype, flxzno3=flxzno3, conczno3=conczno3, &
                                 & dum_swiflux_NO3=dum_swiflux_NO3, dum_POC_conc_swi=dum_POC_conc_swi, &
                                 & dum_RCM_approx=dum_RCM_approx, dum_z_vector=dum_z_vector, &
                                 & dum_diag_profile=dum_diag_profile)

dum_swiflux_NO3 = dum_swiflux_NO3 - por*w*(dum_swiconc_NO3-conczinf)

END SUBROUTINE OMENSED_calculate_zNO3

!========================================================================

SUBROUTINE OMENSED_calculate_zNO3_calcbc(zNO3, bctype, flxzno3, conczno3, dum_swiflux_NO3, &
                                       & dum_POC_conc_swi, dum_RCM_approx, dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zNO3_calcbc* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED_calculate_zNO3 / FUN_OMENSED_zNO3
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: zNO3, flxzno3, conczNO3, dum_swiflux_NO3
INTEGER :: bctype
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAL(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
!
INTEGER :: ltype1, ltype2
REAL(8) :: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
REAL(8) :: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
REAL(8) :: e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
REAL(8) :: e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
REAL(8) :: e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
REAL(8) :: zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
REAL(8) :: e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00
REAL(8) :: e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
REAL(8) :: bctype1_A2, bctype1_B2, bctype2_A2, bctype2_B2
REAL(8) :: rNO3_A2, rNO3_B2, rNO3_A1, rNO3_B1
REAL(8) :: FNH4
INTEGER :: j

! Calculate trial solution for given zno3, matching boundary conditions from layer-by-layer solutions
! Preparation: for each layer, sort out solution - matching across bioturbation boundary (if necessary)

! layer 1: 0 < z < zox, nitrification
CALL sub_prepfg_l12(reac=gammaNH4*NC1, ktemp=0.0D0, zU=0.0D0, zL=zox, D1=DN1, D2=DN2, ls_a=ls_a1, ls_b=ls_b1, &
                  & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ltype=ltype1, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 2: zox < z < zno3, denitrification
CALL sub_prepfg_l12(reac=-NO3CR, ktemp=0.0D0, zU=zox, zL=zNO3, D1=DN1, D2=DN2, ls_a=ls_a2, ls_b=ls_b2, &
                  & ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, ltype=ltype2, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 2 zno3
CALL sub_calcfg_l12(z=zno3, reac=-NO3CR, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ls_D1=DN1, ls_D2=DN2, ltype=ltype2, e=e2_zno3, dedz=dedz2_zno3, f=f2_zno3, dfdz=dfdz2_zno3, &
                  & g=g2_zno3, dgdz=dgdz2_zno3, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                  
! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 -> NO3 source)
! basis functions at bottom of layer 1
CALL sub_calcfg_l12(z=zox, reac=gammaNH4*NC1, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ls_D1=DN1, ls_D2=DN2, ltype=ltype1, e=e1_zox, dedz=dedz1_zox, f=f1_zox, dfdz=dfdz1_zox, &
                  & g=g1_zox, dgdz=dgdz1_zox, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)


! basis functions at top of layer 2
CALL sub_calcfg_l12(z=zox, reac=-NO3CR, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ls_D1=DN1, ls_D2=DN2, ltype=ltype2, e=e2_zox, dedz=dedz2_zox, f=f2_zox, dfdz=dfdz2_zox, &
                  & g=g2_zox, dgdz=dgdz2_zox, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)


! flux of NH4 to zox  TODO NH4 production by denitrification?
FNH4 = FUN_OMENSED_calcReac(zU=zno3, zL=zinf, reac=NC1/(1.0D0+KNH4), &
                          & dum_RCM_approx=dum_RCM_approx, dum_POC_conc_swi=dum_POC_conc_swi)   ! MULTIPLY BY 1/POR ????

! match solutions at zox - continuous concentration, flux discontinuity from H2S ox
IF (zox.LE.zbio) THEN
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, -r_zxf*gammaNH4*FNH4/DN1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
ELSE
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, -r_zxf*gammaNH4*FNH4/DN2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
END IF

! Solution at swi, top of layer 1
CALL sub_calcfg_l12(z=0.0D0, reac=gammaNH4*NC1, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ls_D1=DN1, ls_D2=DN2, ltype=ltype1, e=e1_00, dedz=dedz1_00, f=f1_00, dfdz=dfdz1_00, &
                  & g=g1_00, dgdz=dgdz1_00, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)                  

! transform to use coeffs from l2
CALL sub_xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

! Solve for ANO3, BNO3 given boundary conditions (expressed in terms of transformed basis fns, layer 2 A, B)

! Case 1 zero concentration at zno3
! ANO3*e2_zno3   +  BNO3*f2_zno3  + g2_zno3 = 0;
! ANO3*e1_0     +   BNO3*f1_0     + g1_0  = swi.dum_swiconc_NO3;

! | e2_zno3 f2_zno3 |  |ANO3|   = | -g2_zno3       |
! | e1_0     f1_0   |  |BNO3|     | swi.dum_swiconc_NO3 - g1_0 |


CALL sub_solve2eqn(e2_zno3, f2_zno3, e1_0, f1_0, -g2_zno3, dum_swiconc_NO3 - g1_0, bctype1_A2, bctype1_B2)

! Case  2 zero flux at zno3
! ANO3*de2dz_zno3   +  BNO3*dfdz2_zno3  + dgdz2_zno3 = 0;
! ANO3*e1_0         +   BNO3*f1_0       + g1_0       = swi.dum_swiconc_NO3;

CALL sub_solve2eqn(dedz2_zno3, dfdz2_zno3, e1_0, f1_0, -dgdz2_zno3, dum_swiconc_NO3 - g1_0, bctype2_A2, bctype2_B2)

! Choose type of solution requested
IF (bctype.EQ.1) THEN
   rNO3_A2 = bctype1_A2
   rNO3_B2 = bctype1_B2
ELSE
   rNO3_A2 = bctype2_A2
   rNO3_B2 = bctype2_B2
END IF

! calculate flux at zno3
IF (zno3.LE.zbio) THEN
   flxzno3 = DN1*(rNO3_A2*dedz2_zno3+rNO3_B2*dfdz2_zno3 + dgdz2_zno3)      ! includes 1/por ie flux per (cm^2 pore area)
ELSE
   flxzno3 = DN2*(rNO3_A2*dedz2_zno3+rNO3_B2*dfdz2_zno3 + dgdz2_zno3)      ! includes 1/por ie flux per (cm^2 pore area)
END IF

conczno3 = rNO3_A2*e2_zno3+rNO3_B2*f2_zno3 + g2_zno3
! flux at swi - DO include por so this is per cm^2 water column area
dum_swiflux_NO3 = por*(DN1*(rNO3_A2*dedz1_0+rNO3_B2*dfdz1_0 + dgdz1_0)) ! - w*dum_swiconc_NO3)              ! NB: use A2, B2 as these are _xformed_ layer 1 basis functions

! save coeffs for layer 1
rNO3_A1 = zox_a*rNO3_A2 + zox_b*rNO3_B2 + zox_e
rNO3_B1 = zox_c*rNO3_A2 + zox_d*rNO3_B2 + zox_f

! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_NO3) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zox) THEN
       
       ! layer 1: 0 < z < zox, nitrification
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=gammaNH4*NC1, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, &
                         & ls_e=ls_e1, ls_f=ls_f1, ls_D1=DN1, ls_D2=DN2, ltype=ltype1, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)                  
   
       dum_diag_profile(j,iarr_OS_NO3) = rNO3_A1*e_zx+rNO3_B1*f_zx + g_zx
   END IF 
   
   IF (dum_z_vector(j).GT.zox .AND. dum_z_vector(j).LE.zno3) THEN
   
       ! layer 2: zox < z < zno3, denitrification
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=-NO3CR, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                         & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DN1, ls_D2=DN2, ltype=ltype2, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
   
       dum_diag_profile(j,iarr_OS_NO3) = rNO3_A2*e_zx+rNO3_B2*f_zx + g_zx
       
   END IF
   
END DO

END SUBROUTINE OMENSED_calculate_zNO3_calcbc

!========================================================================

FUNCTION FUN_OMENSED_zNO3(z)
!************************************************************************
!
! *FUN_OMENSED_zNO3* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED_calculate_zNO3 
!
!************************************************************************
        
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: z
!
!*Local variables
!   
REAL(8) :: FUN_OMENSED_zNO3, flxzno3, conczno3, flxswi

CALL OMENSED_calculate_zNO3_calcbc(zNO3=z, bctype=1, flxzno3=flxzno3, conczno3=conczno3, &
                                 & dum_swiflux_NO3=flxswi, dum_POC_conc_swi=OS_POC_conc_swi, &
                                 & dum_RCM_approx=OS_RCM_array, dum_z_vector=OS_z_vector, &
                                 & dum_diag_profile=OS_diag_profile)

FUN_OMENSED_zNO3 = -flxzno3

END FUNCTION FUN_OMENSED_zNO3

!========================================================================

SUBROUTINE OMENSED_calculate_zSO4(dum_swiconc_SO4, dum_swiflux_SO4, dum_POC_conc_swi, dum_RCM_approx, &
                             & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zSO4* calculate sulfate penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
        
IMPLICIT NONE
!
!*Arguments
! 
REAL(8), INTENT(IN) :: dum_swiconc_SO4                ! SO4 concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_SO4         ! SO4 flux: TODO check! (+) sediment -> bottom waters
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAL(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
! 
REAL(8) :: flxzso4, conczso4, loc_conczinf, zL, tol
INTEGER :: bctype

print*, '---------------------- START zSO4 ------------------------------- '

! Iteratively solve for zso4

! try zero flux at zinf and see if we have any SO4 left, also
! calculate [SO4] at zinf for advective loss
bctype = 2
CALL OMENSED_calculate_zSO4_calcbc(zso4=zinf, bctype=bctype, flxzso4=flxzso4, conczso4=loc_conczinf, &
                                 & flxswi=dum_swiflux_SO4, dum_POC_conc_swi=dum_POC_conc_swi, &
                                 & dum_RCM_approx=dum_RCM_approx, dum_z_vector=dum_z_vector, &
                                 & dum_diag_profile=dum_diag_profile)

IF (zno3.EQ.zinf) THEN
   zso4 = zinf
   bctype = 2
ELSE
   IF (loc_conczinf.GE.0) THEN
      zso4 = zinf
      bctype = 2
   ELSE
      bctype = 1
      loc_conczinf = 0.0
      zL=1D-10
      tol=1D-16
      zso4 = FUN_zbrent(FUN_OMENSED_zSO4, max(zno3,zL), zinf, tol)
   END IF
END IF !(zno3 == zinf)

CALL OMENSED_calculate_zSO4_calcbc(zso4=zso4, bctype=bctype, flxzso4=flxzso4, conczso4=conczso4, &
                                 & flxswi=dum_swiflux_SO4, dum_POC_conc_swi=dum_POC_conc_swi, &
                                 & dum_RCM_approx=dum_RCM_approx, dum_z_vector=dum_z_vector, &
                                 & dum_diag_profile=dum_diag_profile)
                                       
dum_swiflux_SO4 = dum_swiflux_SO4 - por*w*(dum_swiconc_SO4 - loc_conczinf)

IF (ABS(dum_swiflux_SO4).LE.const_real_nullsmall) THEN
   dum_swiflux_SO4 = 0.0
END IF

END SUBROUTINE OMENSED_calculate_zSO4

!========================================================================

SUBROUTINE OMENSED_calculate_zSO4_calcbc(zso4, bctype, flxzso4, conczso4, flxswi, &
                                       & dum_POC_conc_swi, dum_RCM_approx, &
                                       & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zSO4_calcbc* calculate sulfate penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED_calculate_zSO4
!
!************************************************************************    
        
IMPLICIT NONE
!
!*Arguments
!        
REAL(8), INTENT(IN) :: zso4
INTEGER, INTENT(IN) :: bctype
REAL(8), INTENT(INOUT) :: flxzso4, conczso4, flxswi
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAL(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
! 
INTEGER :: ltype1, ltype2, ltype3
REAL(8) :: reac_so4
REAL(8) :: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
REAL(8) :: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
REAL(8) :: ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
REAL(8) :: e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
REAL(8) :: e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
REAL(8) :: e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
REAL(8) :: zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
REAL(8) :: e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
REAL(8) :: e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0
REAL(8) :: e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
REAL(8) :: zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
REAL(8) :: e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
REAL(8) :: e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00
REAL(8) :: bctype1_A3, bctype1_B3, bctype2_A3, bctype2_B3
REAL(8) :: rSO4_A1, rSO4_B1, rSO4_A2, rSO4_B2, rSO4_A3, rSO4_B3
REAL(8) :: FH2S
INTEGER :: j

reac_so4=-SO4C

! Calculate trial solution for given zso4, matching boundary conditions from layer-by-layer solutions

! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
! layer 1: 0 < z < zox, passive diffn
!      ls =      sub_prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)

CALL sub_prepfg_l12(reac=0.0D0, ktemp=0.0D0, zU=0.0D0, zL=zox, D1=DSO41, D2=DSO42, ls_a=ls_a1, ls_b=ls_b1, &
                  & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ltype=ltype1, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 2: zox < z < zno3, passive diffn
CALL sub_prepfg_l12(reac=0.0D0, ktemp=0.0D0, zU=zox, zL=zno3, D1=DSO41, D2=DSO42, ls_a=ls_a2, ls_b=ls_b2, &
                  & ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, ltype=ltype2, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 3: zno3 < z < zso4, SO4 consumption by OM oxidation
CALL sub_prepfg_l12(reac=reac_so4, ktemp=0.0D0, zU=zno3, zL=zso4, D1=DSO41, D2=DSO42, ls_a=ls_a3, ls_b=ls_b3, &
                  & ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, ltype=ltype3, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 3 zso4
CALL sub_calcfg_l12(z=zso4, reac=reac_so4, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, &
                  & ls_f=ls_f3, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype3, e=e3_zso4, dedz=dedz3_zso4, f=f3_zso4, dfdz=dfdz3_zso4, &
                  & g=g3_zso4, dgdz=dgdz3_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                  
! Match at zno3, layer 2 - layer 3 (continuity and flux)
! basis functions at bottom of layer 2
CALL sub_calcfg_l12(z=zno3, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype2, e=e2_zno3, dedz=dedz2_zno3, f=f2_zno3, dfdz=dfdz2_zno3, &
                  & g=g2_zno3, dgdz=dgdz2_zno3, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 3
CALL sub_calcfg_l12(z=zno3, reac=reac_so4, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, &
                  & ls_f=ls_f3, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype3, e=e3_zno3, dedz=dedz3_zno3, f=f3_zno3, dfdz=dfdz3_zno3, &
                  & g=g3_zno3, dgdz=dgdz3_zno3, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! match solutions at zno3 - continuous concentration and flux
CALL sub_matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                 & 0.0D0, 0.0D0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
! flux of H2S to oxic interface (Source of SO4)
! NB: include methane region as AOM will produce sulphide as well..

FH2S = FUN_OMENSED_calcReac(zU=zno3, zL=zso4, reac=SO4C, &
                          & dum_RCM_approx=dum_RCM_approx, &
                          & dum_POC_conc_swi=dum_POC_conc_swi) & ! MULTIPLY BY 1/POR ????
     & + gammaCH4*FUN_OMENSED_calcReac(zU=zso4, zL=zinf, reac=MC, &
                          & dum_RCM_approx=dum_RCM_approx, &
                          & dum_POC_conc_swi=dum_POC_conc_swi)
                                                           
! basis functions at bottom of layer 1
CALL sub_calcfg_l12(z=zox, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype1, e=e1_zox, dedz=dedz1_zox, f=f1_zox, dfdz=dfdz1_zox, &
                  & g=g1_zox, dgdz=dgdz1_zox, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! basis functions at top of layer 2
CALL sub_calcfg_l12(z=zox, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype2, e=e2_zox0, dedz=dedz2_zox0, f=f2_zox0, dfdz=dfdz2_zox0, &
                  & g=g2_zox0, dgdz=dgdz2_zox0, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! transform to use coeffs from l3
CALL sub_xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
                 & zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
                 & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox )

! match solutions at zox - continuous concentration, flux discontinuity from H2S ox
IF (zox.LE.zbio) THEN
    CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                     & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                     & 0.0D0, -r_zxf*gammaH2S*FH2S/DSO41, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
ELSE
    CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                     & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                     & 0.0D0, -r_zxf*gammaH2S*FH2S/DSO42, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
END IF

! Solution at swi, top of layer 1
CALL sub_calcfg_l12(z=0.0D0, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype1, e=e1_00, dedz=dedz1_00, f=f1_00, dfdz=dfdz1_00, &
                  & g=g1_00, dgdz=dgdz1_00, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! transform to use coeffs from l3
CALL sub_xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

! Find solutions for two possible types of lower bc
!  case 1  zero concentration at zso4
! Solve for ASO4, BSO4 given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
! ASO4*e3_zso4   +  BSO4*f3_zso4  + g3_zso4 = 0;
! ASO4*e1_0     +   BSO4*f1_0     + g1_0  = swi.dum_swiconc_SO4

! | e3_zso4 f3_zso4 |  |ASO4|   = | -g3_zso4       |
! | e1_0     f1_0   |  |BSO4|     | swi.dum_swiconc_SO4 - g1_0 |
CALL sub_solve2eqn(e3_zso4, f3_zso4, e1_0, f1_0, -g3_zso4, dum_swiconc_SO4 - g1_0, bctype1_A3, bctype1_B3)

! case  2 zero flux at zso4
! ASO4*de3dz_zso4   +  BSO4*dfdz3_zso4  + dgdz3_zso4 = 0;
! ASO4*e1_0         +   BSO4*f1_0       + g1_0       = swi.dum_swiconc_SO4;
CALL sub_solve2eqn(dedz3_zso4, dfdz3_zso4, e1_0, f1_0, -dgdz3_zso4, dum_swiconc_SO4 - g1_0, bctype2_A3, bctype2_B3)

! Choose type of solution requested
IF (bctype.EQ.1) THEN
   rSO4_A3=bctype1_A3
   rSO4_B3=bctype1_B3
ELSE
   rSO4_A3=bctype2_A3
   rSO4_B3=bctype2_B3
END IF

! calculate conc and flux at zso4
conczso4 = rSO4_A3*e3_zso4+rSO4_B3*f3_zso4 + g3_zso4
IF (zso4.LE.zbio) THEN
   flxzso4 = DSO41*(rSO4_A3*dedz3_zso4+rSO4_B3*dfdz3_zso4 + dgdz3_zso4)        ! includes 1/por ie flux per (cm^2 pore area)
ELSE
   flxzso4 = DSO42*(rSO4_A3*dedz3_zso4+rSO4_B3*dfdz3_zso4 + dgdz3_zso4)
END IF

! flux at swi - DO include por so this is per cm^2 water column area
flxswi = por*(DSO41*(rSO4_A3*dedz1_0+rSO4_B3*dfdz1_0 + dgdz1_0)) ! - w*dum_swiconc_SO4)   ! NB: use A3, B3 as these are _xformed_ layer 1 basis functions

! save coeffs for layers 2 and 1
rSO4_A2 = zno3_a*rSO4_A3 + zno3_b*rSO4_B3 + zno3_e
rSO4_B2 = zno3_c*rSO4_A3 + zno3_d*rSO4_B3 + zno3_f
rSO4_A1 = zox_a*rSO4_A3 + zox_b*rSO4_B3 + zox_e
rSO4_B1 = zox_c*rSO4_A3 + zox_d*rSO4_B3 + zox_f

! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_SO4) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zox) THEN
       
       ! layer 1: 0 < z < zox
       
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=0.0D0, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, &
                         & ls_e=ls_e1, ls_f=ls_f1, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype1, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
   
       dum_diag_profile(j,iarr_OS_SO4) = rSO4_A1*e_zx+rSO4_B1*f_zx + g_zx
   
   END IF
   
   IF (dum_z_vector(j).GT.zox .AND. dum_z_vector(j).LE.zno3) THEN
       
       ! layer 2: zox < z < zno3, denitrification
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                         & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype2, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_SO4) = rSO4_A2*e_zx+rSO4_B2*f_zx + g_zx
   
   END IF
   
   IF (dum_z_vector(j).GT.zno3 .AND. dum_z_vector(j).LE.zso4) THEN
       
       ! layer 3: zno3 < z < zso4, SO4 consumption by OM oxidation
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac_so4, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, &
                         & ls_e=ls_e3, ls_f=ls_f3, ls_D1=DSO41, ls_D2=DSO42, ltype=ltype3, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_SO4) = rSO4_A3*e_zx+rSO4_B3*f_zx + g_zx
   
   END IF 
      
END DO

END SUBROUTINE OMENSED_calculate_zSO4_calcbc

!========================================================================

FUNCTION FUN_OMENSED_calcFSO4(z)
!************************************************************************
!
! *FUN_OMENSED_calcFSO4* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_zSO4
!
!************************************************************************
         
IMPLICIT NONE
!
!*Arguments
!        
REAL(8) :: z
!
!*Local variables
! 
REAL(8) :: FUN_OMENSED_calcFSO4, tmpreac

! Calculate SO4 consumption below zso4, by organic matter and indirectly via methane oxidation
tmpreac    = MC*gammaCH4

FUN_OMENSED_calcFSO4 = FUN_OMENSED_calcReac(zU=z, zL=zinf, reac=tmpreac, dum_RCM_approx=OS_RCM_array, &
                                          & dum_POC_conc_swi=OS_POC_conc_swi)

END FUNCTION FUN_OMENSED_calcFSO4

!========================================================================

FUNCTION FUN_OMENSED_zSO4(z)
!************************************************************************
!
! *FUN_OMENSED_zSO4* calculate oxygen penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED_calculate_zSO4 
!
!************************************************************************
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8) :: z
!
!*Local variables
! 
REAL(8) :: FUN_OMENSED_zSO4, flxzso4, conczso4, flxswi

CALL OMENSED_calculate_zSO4_calcbc(zso4=z, bctype=1, flxzso4=flxzso4, conczso4=conczso4, flxswi=flxswi, &
                                 & dum_POC_conc_swi=OS_POC_conc_swi, dum_RCM_approx=OS_RCM_array, &
                                 & dum_z_vector=OS_z_vector, dum_diag_profile=OS_diag_profile)
                                
FUN_OMENSED_zSO4 = -flxzso4 - FUN_OMENSED_calcFSO4(z)

END FUNCTION FUN_OMENSED_zSO4

!========================================================================

SUBROUTINE OMENSED_calculate_zNH4(dum_swiconc_NH4, dum_swiflux_NH4, dum_POC_conc_swi, dum_RCM_approx, &
                                & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zNH4* calculate ammonium penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: dum_swiconc_NH4                ! NH4 concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_NH4         ! NH4 flux: TODO check! (+) sediment -> bottom waters
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAL(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
! 
REAL(8) :: loc_conczinf
INTEGER :: ltype1, ltype2, ltype3
REAL(8) :: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
REAL(8) :: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
REAL(8) :: ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
REAL(8) :: e3_zinf, dedz3_zinf, f3_zinf, dfdz3_zinf, g3_zinf, dgdz3_zinf
REAL(8) :: e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
REAL(8) :: e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
REAL(8) :: zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
REAL(8) :: e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
REAL(8) :: e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
REAL(8) :: e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
REAL(8) :: zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
REAL(8) :: e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
REAL(8) :: e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
REAL(8) :: rNH4_A3, rNH4_B3
REAL(8) :: rNH4_A2, rNH4_B2
REAL(8) :: rNH4_A1, rNH4_B1
REAL(8) :: FNH4
INTEGER :: j

print*, '---------------------- START zNH4 ------------------------------- '

! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
! layer 1: 0 < z < zox, NH4 prod (remaining after oxidation)
!      ls =      sub_prepfg_l12( bsd, swi, r, reac1,     reac2,     ktemp, zU, zL, D1,        D2)
CALL sub_prepfg_l12(reac=(1.0D0-gammaNH4)*NC1/(1.0D0+KNH4), ktemp=0.0D0, zU=0.0D0, zL=zox, D1=DNH41, &
                  & D2=DNH42, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ltype=ltype1, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 2: zox < z < zno3, passive diffn TODO NH4 from denitrification?
CALL sub_prepfg_l12(reac=0.0D0, ktemp=0.0D0, zU=zox, zL=zno3, D1=DNH41, &
                  & D2=DNH42, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ltype=ltype2, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 3: zno3 < z < zinf, NH4 production
CALL sub_prepfg_l12(reac=NC1/(1.0D0+KNH4), ktemp=0.0D0, zU=zno3, zL=zinf, D1=DNH41, &
                  & D2=DNH42, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, &
                  & ls_f=ls_f3, ltype=ltype3, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 3 zinf
CALL sub_calcfg_l12(z=zinf, reac=NC1/(1.0D0+KNH4), ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, &
                  & ls_e=ls_e3, ls_f=ls_f3, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype3, &
                  & e=e3_zinf, dedz=dedz3_zinf, f=f3_zinf, dfdz=dfdz3_zinf, g=g3_zinf, dgdz=dgdz3_zinf, &
                  & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! Match at zno3, layer 2 - layer 3 (continuity and flux)
! basis functions at bottom of layer 2
CALL sub_calcfg_l12(z=zno3, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                  & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype2, &
                  & e=e2_zno3, dedz=dedz2_zno3, f=f2_zno3, dfdz=dfdz2_zno3, g=g2_zno3, dgdz=dgdz2_zno3, &
                  & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 3
CALL sub_calcfg_l12(z=zno3, reac=NC1/(1.0D0+KNH4), ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, &
                  & ls_e=ls_e3, ls_f=ls_f3, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype3, &
                  & e=e3_zno3, dedz=dedz3_zno3, f=f3_zno3, dfdz=dfdz3_zno3, g=g3_zno3, dgdz=dgdz3_zno3, &
                  & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! match solutions at zno3 - continuous concentration and flux
CALL sub_matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                 & 0.0D0, 0.0D0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from NH4 sink)
! flux of NH4 to oxic interface  TODO NH4 prod by denitrification?
FNH4 = FUN_OMENSED_calcReac(zU=zno3, zL=zinf, reac=NC1/(1.0D0+KNH4), &
                          & dum_RCM_approx=dum_RCM_approx, dum_POC_conc_swi=dum_POC_conc_swi)   ! MULTIPLY BY 1/POR ????
                                         
! basis functions at bottom of layer 1
CALL sub_calcfg_l12(z=zox, reac=(1.0D0-gammaNH4)*NC1/(1.0D0+KNH4), ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, &
                  & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype1, &
                  & e=e1_zox, dedz=dedz1_zox, f=f1_zox, dfdz=dfdz1_zox, g=g1_zox, dgdz=dgdz1_zox, &
                  & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! basis functions at top of layer 2
CALL sub_calcfg_l12(z=zox, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, &
                  & ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype2, &
                  & e=e2_zox0, dedz=dedz2_zox0, f=f2_zox0, dfdz=dfdz2_zox0, g=g2_zox0, dgdz=dgdz2_zox0, &
                  & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! transform to use coeffs from l3
CALL sub_xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
                 & zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
                 & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

! match solutions at zox - continuous concentration, flux discontinuity from NH4 ox
IF (zox.LE.zbio) THEN
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, r_zxf*gammaNH4*FNH4/DNH41, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
ELSE
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, r_zxf*gammaNH4*FNH4/DNH42, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
END IF

! Solution at swi, top of layer 1
CALL sub_calcfg_l12(z=0.0D0, reac=(1.0D0-gammaNH4)*NC1/(1.0D0+KNH4), ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, &
                  & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype1, &
                  & e=e1_00, dedz=dedz1_00, f=f1_00, dfdz=dfdz1_00, g=g1_00, dgdz=dgdz1_00, &
                  & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! transform to use coeffs from l3
CALL sub_xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
                 & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

! Solve for ANH4, BNH4 given boundary conditions (expressed in terms of transformed basis fns, layer 3 A, B)
! ANH4*dedz3_zinf   +  BNH4*dfdz3_zinf  + dgdz3_zinf = 0;
! ANH4*e1_0     +   BNH4*f1_0     + g1_0  = swi.;
! | dedz3_zinf dfdz3_zinf |  |ANH4|   = | -dgdz3_zinf       |
! | e1_0     f1_0         |  |BNH4|     | swi. - g1_0 |

CALL sub_solve2eqn(dedz3_zinf, dfdz3_zinf, e1_0, f1_0, -dgdz3_zinf,  - g1_0, rNH4_A3, rNH4_B3)

! calculate concentration at zinf
loc_conczinf = rNH4_A3*e3_zinf+rNH4_B3*f3_zinf + g3_zinf

! flux at swi - DO include por so this is per cm^2 water column area
dum_swiflux_NH4 = por*(DNH41*(rNH4_A3*dedz1_0+rNH4_B3*dfdz1_0 + dgdz1_0) - w*(dum_swiconc_NH4 - loc_conczinf))   ! NB: use A3, B3 as these are _xformed_ layer 1 basis functions

IF (ABS(dum_swiflux_NH4).LE.const_real_nullsmall) THEN
   dum_swiflux_NH4 = 0.0
END IF

! save coeffs for layers 2 and 1
rNH4_A2 = zno3_a*rNH4_A3 + zno3_b*rNH4_B3 + zno3_e
rNH4_B2 = zno3_c*rNH4_A3 + zno3_d*rNH4_B3 + zno3_f
rNH4_A1 = zox_a*rNH4_A3 + zox_b*rNH4_B3 + zox_e
rNH4_B1 = zox_c*rNH4_A3 + zox_d*rNH4_B3 + zox_f

! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_NH4) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zox) THEN
       
       ! layer 1: 0 < z < zox
       
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=(1.0D0-gammaNH4)*NC1/(1.0D0+KNH4), ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, &
                         & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype1, &
                         & e=e_zx, dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
   
       dum_diag_profile(j,iarr_OS_NH4) = rNH4_A1*e_zx+rNH4_B1*f_zx + g_zx
       
   END IF
   
   IF (dum_z_vector(j).GT.zox .AND. dum_z_vector(j).LE.zno3) THEN
       
       ! layer 2: zox < z < zno3, denitrification
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                         & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype2, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_NH4) = rNH4_A2*e_zx+rNH4_B2*f_zx + g_zx
   
   END IF
   
   IF (dum_z_vector(j).GT.zno3 ) THEN
       
       ! layer 3: zno3 < z < zinf, NH4 production
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=NC1/(1.0D0+KNH4), ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, &
                         & ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, ls_D1=DNH41, ls_D2=DNH42, ltype=ltype3, e=e_zx, &
                         & dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_NH4) = rNH4_A3*e_zx+rNH4_B3*f_zx + g_zx
   
   END IF 
      
END DO

END SUBROUTINE OMENSED_calculate_zNH4

!========================================================================

SUBROUTINE OMENSED_calculate_zH2S(dum_swiconc_H2S, dum_swiflux_H2S, &
                                & dum_POC_conc_swi, dum_RCM_approx, &
                                & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zH2S* calculate ammonium penetration depth 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: dum_swiconc_H2S                ! SO4 concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_H2S         ! SO4 flux: TODO check! (+) sediment -> bottom waters
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
! 
REAL(8) :: loc_conczinf
REAL(8) :: reac_h2s                ! reactive terms: OM degradation
INTEGER :: ltype1, ltype2, ltype3, ltype4
REAL(8) :: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
REAL(8) :: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
REAL(8) :: ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
REAL(8) :: ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4
REAL(8) :: e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf
REAL(8) :: e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
REAL(8) :: e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4
REAL(8) :: zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f
REAL(8) :: e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
REAL(8) :: e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30
REAL(8) :: e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
REAL(8) :: zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
REAL(8) :: e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
REAL(8) :: e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
REAL(8) :: e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
REAL(8) :: zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
REAL(8) :: e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
REAL(8) :: e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
REAL(8) :: rH2S_A4, rH2S_B4
REAL(8) :: rH2S_A3, rH2S_B3
REAL(8) :: rH2S_A2, rH2S_B2
REAL(8) :: rH2S_A1, rH2S_B1
REAL(8) :: zso4FH2S, zoxFH2S
INTEGER :: j
    
reac_h2s=SO4C

print*, '---------------------- START zH2S ------------------------------- '
                      
! Calculate H2S
! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
! layer 1: 0 < z < zox, passive diffn
CALL sub_prepfg_l12(reac=0.0D0, ktemp=0.0D0, zU=0.0D0, zL=zox, D1=DH2S1, &
                  & D2=DH2S2, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ltype=ltype1, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 2: zox < z < zno3, passive diffn
CALL sub_prepfg_l12(reac=0.0D0, ktemp=0.0D0, zU=zox, zL=zno3, D1=DH2S1, &
                  & D2=DH2S2, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ltype=ltype2, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 3: zno3 < z < zso4, H2S consumption by oxidation
CALL sub_prepfg_l12(reac=reac_h2s, ktemp=0.0D0, zU=zno3, zL=zso4, D1=DH2S1, &
                  & D2=DH2S2, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, &
                  & ls_f=ls_f3, ltype=ltype3, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 4: zso4 < z < zinf, passive diffn
CALL sub_prepfg_l12(reac=0.0D0, ktemp=0.0D0, zU=zso4, zL=zinf, D1=DH2S1, &
                  & D2=DH2S2, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, ls_e=ls_e4, &
                  & ls_f=ls_f4, ltype=ltype4, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 4 zinf
CALL sub_calcfg_l12(z=zinf, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, ls_e=ls_e4, ls_f=ls_f4, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype4, e=e4_zinf, dedz=dedz4_zinf, f=f4_zinf, dfdz=dfdz4_zinf, &
                  & g=g4_zinf, dgdz=dgdz4_zinf, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
! basis functions at bottom of layer 3
CALL sub_calcfg_l12(z=zso4, reac=reac_h2s, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype3, e=e3_zso4, dedz=dedz3_zso4, f=f3_zso4, dfdz=dfdz3_zso4, &
                  & g=g3_zso4, dgdz=dgdz3_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 4
CALL sub_calcfg_l12(z=zso4, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, ls_e=ls_e4, ls_f=ls_f4, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype4, e=e4_zso4, dedz=dedz4_zso4, f=f4_zso4, dfdz=dfdz4_zso4, &
                  & g=g4_zso4, dgdz=dgdz4_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! flux of H2S produced by AOM interface (Source of H2S)
zso4FH2S = FUN_OMENSED_calcReac(zU=zso4, zL=zinf, reac=MC, dum_RCM_approx=dum_RCM_approx, &
                              & dum_POC_conc_swi=dum_POC_conc_swi) ! MULTIPLY BY 1/POR ????
                                          
! match solutions at zso4 - continuous concentration and flux
CALL sub_matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, &
                 & e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, &
                 & 0.0D0, -gammaCH4*zso4FH2S/DH2S2, zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f)

! Match at zno3, layer 2 - layer 3 (continuity and flux)
! basis functions at bottom of layer 2
CALL sub_calcfg_l12(z=zno3, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype2, e=e2_zno3, dedz=dedz2_zno3, f=f2_zno3, dfdz=dfdz2_zno3, &
                  & g=g2_zno3, dgdz=dgdz2_zno3, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 3
CALL sub_calcfg_l12(z=zno3, reac=reac_h2s, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype3, e=e3_zno30, dedz=dedz3_zno30, f=f3_zno30, dfdz=dfdz3_zno30, &
                  & g=g3_zno30, dgdz=dgdz3_zno30, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... transformed to use coeffs from l4
CALL sub_xformsoln(e3_zno30, f3_zno30, g3_zno30, dedz3_zno30, dfdz3_zno30, dgdz3_zno30, &
                 & zso4_a , zso4_b , zso4_c , zso4_d , zso4_e ,zso4_f, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3)

! match solutions at zno3 - continuous concentration and flux
CALL sub_matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                 & 0.0D0, 0.0D0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from H2S source)
! flux of H2S to oxic interface (from all sources of H2S below)
! NB: include methane region as AOM will produce sulphide as well..
zoxFH2S = FUN_OMENSED_calcReac(zU=zno3, zL=zso4, reac=SO4C, &
                             & dum_RCM_approx=dum_RCM_approx, &
                             & dum_POC_conc_swi=dum_POC_conc_swi) + &
        & FUN_OMENSED_calcReac(zU=zso4, zL=zinf, reac=MC, &
                             & dum_RCM_approx=dum_RCM_approx, &
                             & dum_POC_conc_swi=dum_POC_conc_swi)

! basis functions at bottom of layer 1
CALL sub_calcfg_l12(z=zox, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype1, e=e1_zox, dedz=dedz1_zox, f=f1_zox, dfdz=dfdz1_zox, &
                  & g=g1_zox, dgdz=dgdz1_zox, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! basis functions at top of layer 2
CALL sub_calcfg_l12(z=zox, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype2, e=e2_zox0, dedz=dedz2_zox0, f=f2_zox0, dfdz=dfdz2_zox0, &
                  & g=g2_zox0, dgdz=dgdz2_zox0, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

!   transform to use coeffs from l4
CALL sub_xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
                 & zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
                 & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

! match solutions at zox - continuous concentration, flux discontinuity from H2S ox
IF (zox.LE.zbio) THEN
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, r_zxf*gammaH2S*zoxFH2S/DH2S1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
ELSE
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, r_zxf*gammaH2S*zoxFH2S/DH2S2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
END IF

! Solution at swi, top of layer 1
CALL sub_calcfg_l12(z=0.0D0, reac=0.0D0, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, &
                  & ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype1, e=e1_00, dedz=dedz1_00, f=f1_00, dfdz=dfdz1_00, &
                  & g=g1_00, dgdz=dgdz1_00, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)


! transform to use coeffs from l4
CALL sub_xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
                 & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

! Solve for AH2S, BH2S given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
!  AH2S*dedz4_zinf   +  BH2S*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
!  AH2S*e1_0     +   BH2S*f1_0     + g1_0  = swi.dum_swiconc_H2S;
!  | dedz4_zinf dfdz4_zinf |  |AH2S|   = | -dgz4_zinf       |
!  | e1_0     f1_0         |  |BH2S|     | swi.dum_swiconc_H2S - g1_0 |
CALL sub_solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, dum_swiconc_H2S - g1_0, rH2S_A4, rH2S_B4)

! calculate concentration at zinf
loc_conczinf = rH2S_A4*e4_zinf+rH2S_B4*f4_zinf + g4_zinf

! flux at swi - DO include por so this is per cm^2 water column area
dum_swiflux_H2S = por*(DH2S1*(rH2S_A4*dedz1_0+rH2S_B4*dfdz1_0 + dgdz1_0) - w*(dum_swiconc_H2S - loc_conczinf))   ! NB: use A4, B4 as these are _xformed_ layer 1 basis functions

IF (ABS(dum_swiflux_H2S).LE.const_real_nullsmall) THEN
    dum_swiflux_H2S = 0.0
END IF

! save coeffs for layers 3, 2 and 1
rH2S_A3 = zso4_a*rH2S_A4 + zso4_b*rH2S_B4 + zso4_e
rH2S_B3 = zso4_c*rH2S_A4 + zso4_d*rH2S_B4 + zso4_f
rH2S_A2 = zno3_a*rH2S_A4 + zno3_b*rH2S_B4 + zno3_e
rH2S_B2 = zno3_c*rH2S_A4 + zno3_d*rH2S_B4 + zno3_f
rH2S_A1 = zox_a*rH2S_A4 + zox_b*rH2S_B4 + zox_e
rH2S_B1 = zox_c*rH2S_A4 + zox_d*rH2S_B4 + zox_f

! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_H2S) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zox) THEN
       
       ! layer 1: 0 < z < zox
       
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=0.0D0, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, &
                         & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype1, &
                         & e=e_zx, dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
   
       dum_diag_profile(j,iarr_OS_H2S) = rH2S_A1*e_zx+rH2S_B1*f_zx + g_zx
       
   END IF
   
   IF (dum_z_vector(j).GT.zox .AND. dum_z_vector(j).LE.zno3) THEN
       
       ! layer 2: zox < z < zno3, denitrification
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=0.0D0, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                         & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype2, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_H2S) = rH2S_A2*e_zx+rH2S_B2*f_zx + g_zx
   
   END IF
   
   IF (dum_z_vector(j).GT.zno3 .AND. dum_z_vector(j).LE.zso4) THEN
       
       ! layer 3: zno3 < z < zso4, SO4 consumption by OM oxidation
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac_h2s, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, &
                         & ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype3, e=e_zx, &
                         & dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_H2S) = rH2S_A3*e_zx+rH2S_B3*f_zx + g_zx
       
   END IF
   
   IF (dum_z_vector(j).GT.zso4) THEN
       
       ! layer 4: zso4 < z < zinf, passive diffn
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=0.0D0, ktemp=0.0D0, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, &
                         & ls_d=ls_d4, ls_e=ls_e4, ls_f=ls_f4, ls_D1=DH2S1, ls_D2=DH2S2, ltype=ltype4, e=e_zx, &
                         & dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_H2S) = rH2S_A4*e_zx+rH2S_B4*f_zx + g_zx
   
   END IF 
      
END DO

END SUBROUTINE OMENSED_calculate_zH2S

!========================================================================

SUBROUTINE OMENSED_calculate_zPO4_M(dum_swiconc_PO4, dum_swiflux_PO4, dum_swiflux_M, &
                                  & dum_POC_conc_swi, dum_RCM_approx, &
                                  & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zPO4_M* calculate phosphate 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!  
REAL(8), INTENT(IN) :: dum_swiconc_PO4            ! PO4 concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_PO4     ! PO4 flux: TODO check! (+) sediment -> bottom waters
REAL(8), INTENT(INOUT) :: dum_swiflux_M       ! PO4 flux: TODO check! (+) sediment -> bottom waters
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8), DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
! 
REAL(8) :: loc_conczinf
REAL(8) :: rPO4_M_A1, rPO4_M_B1, rPO4_M_C1, rPO4_M_D1
REAL(8) :: rPO4_M_A2, rPO4_M_B2, rPO4_M_C2, rPO4_M_D2
INTEGER :: dum_ltype1, dum_ltype2
REAL(8), DIMENSION(4,4) :: dum_mat_C1, dum_mat_C2
REAL(8), DIMENSION(1:4) ::  dum_vec_D1, dum_vec_D2
INTEGER :: loc_i, loc_j
REAL(8) :: reac1_po4_ox, reac2_po4_ox, reac1_po4_anox, reac2_po4_anox                 ! reactive terms: OM degradation
REAL(8) :: e2_zinf_P, dedz2_zinf_P, f2_zinf_P, dfdz2_zinf_P, g2_zinf_P, dgdz2_zinf_P
REAL(8) :: p2_zinf_P, dpdz2_zinf_P, q2_zinf_P, dqdz2_zinf_P
REAL(8) :: e2_zinf_M, dedz2_zinf_M, f2_zinf_M, dfdz2_zinf_M, g2_zinf_M, dgdz2_zinf_M
REAL(8) :: p2_zinf_M, dpdz2_zinf_M, q2_zinf_M, dqdz2_zinf_M
REAL(8) :: e1_zox_P, dedz1_zox_P, f1_zox_P, dfdz1_zox_P, g1_zox_P, dgdz1_zox_P
REAL(8) :: p1_zox_P, dpdz1_zox_P, q1_zox_P, dqdz1_zox_P
REAL(8) :: e1_zox_M, dedz1_zox_M, f1_zox_M, dfdz1_zox_M, g1_zox_M, dgdz1_zox_M
REAL(8) :: p1_zox_M, dpdz1_zox_M, q1_zox_M, dqdz1_zox_M
REAL(8) :: e2_zox_P, dedz2_zox_P, f2_zox_P, dfdz2_zox_P, g2_zox_P, dgdz2_zox_P
REAL(8) :: p2_zox_P, dpdz2_zox_P, q2_zox_P, dqdz2_zox_P
REAL(8) :: e2_zox_M, dedz2_zox_M, f2_zox_M, dfdz2_zox_M, g2_zox_M, dgdz2_zox_M
REAL(8) :: p2_zox_M, dpdz2_zox_M, q2_zox_M, dqdz2_zox_M
REAL(8) :: e1_z0_P, dedz1_z0_P, f1_z0_P, dfdz1_z0_P, g1_z0_P, dgdz1_z0_P
REAL(8) :: p1_z0_P, dpdz1_z0_P, q1_z0_P, dqdz1_z0_P
REAL(8) :: e1_z0_M, dedz1_z0_M, f1_z0_M, dfdz1_z0_M, g1_z0_M, dgdz1_z0_M
REAL(8) :: p1_z0_M, dpdz1_z0_M, q1_z0_M, dqdz1_z0_M
REAL(8) :: e_zx_P, dedz_zx_P, f_zx_P, dfdz_zx_P, g_zx_P, dgdz_zx_P
REAL(8) :: p_zx_P, dpdz_zx_P, q_zx_P, dqdz_zx_P
REAL(8) :: e_zx_M, dedz_zx_M, f_zx_M, dfdz_zx_M, g_zx_M, dgdz_zx_M
REAL(8) :: p_zx_M, dpdz_zx_M, q_zx_M, dqdz_zx_M
REAL(8) :: g_P, dgdz_P
REAL(8) :: g_M, dgdz_M
REAL(8) :: loc_Vb, loc_Fb                                   ! discontinuity constants
REAL(8), DIMENSION(1:4,1:4) :: loc_mat_X_4x4, loc_mat_Y_4x4
REAL(8), DIMENSION(1:4) ::  loc_vec_Z_4
REAL(8), DIMENSION(1:3,1:3) :: loc_mat_X_3x3, loc_mat_Y_3x3
REAL(8), DIMENSION(1:3) ::  loc_vec_Z_3
REAL(8), DIMENSION(4,4) :: loc_mat_C_4x4
REAL(8), DIMENSION(1:4) ::  loc_vec_D_4
REAL(8), DIMENSION(3,3) :: loc_mat_C_3x3
REAL(8), DIMENSION(1:3) ::  loc_vec_D_3
INTEGER :: loc_dim
REAL(8), DIMENSION(1:4) :: loc_EFPQ_P, loc_dEFPQdz_P, loc_EFPQ_M, loc_dEFPQdz_M
REAL(8), DIMENSION(1:4) :: loc_EFPQ_P_t, loc_dEFPQdz_P_t, loc_EFPQ_M_t, loc_dEFPQdz_M_t
REAL(8), DIMENSION(1:4) :: loc_Layer1_IC, loc_Layer2_IC
INTEGER :: j

reac1_po4_ox = 1.0D0/(1.0D0+KPO4_ox)*PC1
reac2_po4_ox = 1.0D0/(1.0D0+KPO4_ox)*PC2
reac1_po4_anox = 1.0D0/(1.0D0+KPO4_anox)*PC1
reac2_po4_anox = 1.0D0/(1.0D0+KPO4_anox)*PC2

loc_Vb = 0.0D0
loc_Fb = 0.0D0

! Initialize loc_mat_C_4x4 & loc_vec_D_4 with zeros as, so I don't need to add them later manually ...
loc_mat_C_4x4 = 0.0D0
loc_vec_D_4 = 0.0D0

print*, '---------------------- START zPO4_M ------------------------------- '

!   Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
!   layer 1: 0 < z < zox, OM degradation (-) Sorption to sediment Fe-oxides (ktemp)
CALL sub_prepfg_l12_PO4_M(reac1=reac1_po4_ox, ktempP=ksPO4/(1.0D0+KPO4_ox), QtempP=PO4s*ksPO4/(1.0D0+KPO4_ox), &
                        & zU=0.0D0, zL=zox, D1P=DPO41/(1.0D0+KPO4_ox), D2P=DPO42/(1.0D0+KPO4_ox), alphaP=0.0D0, &
                        & ktempM=0.0D0, QtempM=0.0D0, D1M=Dbio, D2M=0.0D0, alphaM=(1.0D0/SD)*ksPO4, loc_mat_C=dum_mat_C1, &
                        & loc_vec_D=dum_vec_D1, ltype=dum_ltype1, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                        
!   layer 2: zox < z < zinf,
!   OM degradation (-) authigenic P formation (ktemp) (+) P desorption due to Fe-bound P release upon Fe oxide reduction
CALL sub_prepfg_l12_PO4_M(reac1=reac1_po4_anox, ktempP=kaPO4/(1.0D0+KPO4_anox), QtempP=PO4a*kaPO4/(1.0D0+KPO4_anox), &
                        & zU=zox, zL=zinf, D1P=DPO41/(1.0D0+KPO4_anox), D2P=DPO42/(1.0D0+KPO4_anox), &
                        & alphaP=SD*kmPO4/(1.0D0+KPO4_anox), &
                        & ktempM=kmPO4, QtempM=kmPO4*Minf, D1M=Dbio, D2M=0.0D0, alphaM=0.0D0, loc_mat_C=dum_mat_C2, &
                        & loc_vec_D=dum_vec_D2, &
                        & ltype=dum_ltype2, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 2 zinf
CALL sub_calcfg_l12_PO4_M(z=zinf, reac1P=reac1_po4_anox, dum_ktempP=kaPO4/(1.0D0+KPO4_anox), &
                        & dum_QtempP=PO4a*kaPO4/(1.0D0+KPO4_anox), &
                        & dum_D1P=DPO41/(1.0D0+KPO4_anox), dum_D2P=DPO42/(1+KPO4_anox), dum_alphaP=SD*kmPO4/(1+KPO4_anox), &
                        & dum_mat_C=dum_mat_C2, dum_vec_D=dum_vec_D2, dum_ltype=dum_ltype2, dum_ktempM=kmPO4, &
                        & dum_QtempM=kmPO4*Minf, &
                        & dum_D1M=Dbio, dum_D2M=0.0D0, dum_alphaM=0.0D0, e_P=e2_zinf_P, dedz_P=dedz2_zinf_P, f_P=f2_zinf_P, &
                        & dfdz_P=dfdz2_zinf_P, g_P=g2_zinf_P, dgdz_P=dgdz2_zinf_P, p_P=p2_zinf_P, dpdz_P=dpdz2_zinf_P, &
                        & q_P=q2_zinf_P, &
                        & dqdz_P=dqdz2_zinf_P, e_M=e2_zinf_M, dedz_M=dedz2_zinf_M, f_M=f2_zinf_M, dfdz_M=dfdz2_zinf_M, &
                        & g_M=g2_zinf_M, &
                        & dgdz_M=dgdz2_zinf_M, p_M=p2_zinf_M, dpdz_M=dpdz2_zinf_M, q_M=q2_zinf_M, dqdz_M=dqdz2_zinf_M, &
                        & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)


! Match at zox, layer 1 - layer 2 (continuity and flux)
! basis functions at bottom of layer 1
CALL sub_calcfg_l12_PO4_M(z=zox, reac1P=reac1_po4_ox, dum_ktempP=ksPO4/(1.0D0+KPO4_ox), dum_QtempP=PO4s*ksPO4/(1.0D0+KPO4_ox), &
                        & dum_D1P=DPO41/(1.0D0+KPO4_ox), dum_D2P=DPO42/(1.0D0+KPO4_ox), dum_alphaP=0.0D0, &
                        & dum_mat_C=dum_mat_C1, dum_vec_D=dum_vec_D1, dum_ltype=dum_ltype1, dum_ktempM=0.0D0, dum_QtempM=0.0D0, &
                        & dum_D1M=Dbio, dum_D2M=0.0D0, dum_alphaM=(1.0D0/SD)*ksPO4, e_P=e1_zox_P, dedz_P=dedz1_zox_P, &
                        & f_P=f1_zox_P, &
                        & dfdz_P=dfdz1_zox_P, g_P=g1_zox_P, dgdz_P=dgdz1_zox_P, p_P=p1_zox_P, dpdz_P=dpdz1_zox_P, &
                        & q_P=q1_zox_P, &
                        & dqdz_P=dqdz1_zox_P, e_M=e1_zox_M, dedz_M=dedz1_zox_M, f_M=f1_zox_M, dfdz_M=dfdz1_zox_M, &
                        & g_M=g1_zox_M, &
                        & dgdz_M=dgdz1_zox_M, p_M=p1_zox_M, dpdz_M=dpdz1_zox_M, q_M=q1_zox_M, dqdz_M=dqdz1_zox_M, &
                        & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

!  and top of layer 2
CALL sub_calcfg_l12_PO4_M(z=zox, reac1P=reac1_po4_anox, dum_ktempP=kaPO4/(1.0D0+KPO4_anox), &
                        & dum_QtempP=PO4a*kaPO4/(1.0D0+KPO4_anox), &
                        & dum_D1P=DPO41/(1.0D0+KPO4_anox), dum_D2P=DPO42/(1.0D0+KPO4_anox), &
                        & dum_alphaP=SD*kmPO4/(1.0D0+KPO4_anox), &
                        & dum_mat_C=dum_mat_C2, dum_vec_D=dum_vec_D2, dum_ltype=dum_ltype2, dum_ktempM=kmPO4, &
                        & dum_QtempM=kmPO4*Minf, &
                        & dum_D1M=Dbio, dum_D2M=0.0D0, dum_alphaM=0.0D0, e_P=e2_zox_P, dedz_P=dedz2_zox_P, f_P=f2_zox_P, &
                        & dfdz_P=dfdz2_zox_P, g_P=g2_zox_P, dgdz_P=dgdz2_zox_P, p_P=p2_zox_P, dpdz_P=dpdz2_zox_P, &
                        & q_P=q2_zox_P, &
                        & dqdz_P=dqdz2_zox_P, e_M=e2_zox_M, dedz_M=dedz2_zox_M, f_M=f2_zox_M, dfdz_M=dfdz2_zox_M, &
                        & g_M=g2_zox_M, &
                        & dgdz_M=dgdz2_zox_M, p_M=p2_zox_M, dpdz_M=dpdz2_zox_M, q_M=q2_zox_M, dqdz_M=dqdz2_zox_M, &
                        & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! match solutions at zox - continuous concentration and flux
! organize the data in matrices and let the intrinsic fortran function do the calculation
! DH: Maybe this could be done more efficiently !?
!  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity
!  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux
!  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity
!  |       x16|   | D_l |      |        y16| | D_r|    |z4|    SD M flux only in bioturbated case, otherwise not an independent constraint

IF (zox.LT.zbio) THEN  ! 1. CASE: 4 int const. in each layer
   ! weird FORTRAN matrices makes the transpose necessary
   loc_mat_X_4x4 = transpose(reshape((/ e1_zox_P, f1_zox_P, p1_zox_P, q1_zox_P, &
                          & dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P, dqdz1_zox_P, &
                          & e1_zox_M, f1_zox_M, p1_zox_M, q1_zox_M, &
                          & dedz1_zox_M, dfdz1_zox_M, dpdz1_zox_M, dqdz1_zox_M /), shape(loc_mat_X_4x4)))

   loc_mat_Y_4x4 = transpose(reshape((/ e2_zox_P, f2_zox_P, p2_zox_P, q2_zox_P, &
                          & dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P, dqdz2_zox_P, &
                          & e2_zox_M, f2_zox_M, p2_zox_M, q2_zox_M, &
                          & dedz2_zox_M, dfdz2_zox_M, dpdz2_zox_M, dqdz2_zox_M /), shape(loc_mat_Y_4x4)))

   loc_vec_Z_4 = (/ g2_zox_P-g1_zox_P + loc_Vb, &
                 & dgdz2_zox_P - dgdz1_zox_P + loc_Fb - w*loc_Vb, &
                 & g2_zox_M-g1_zox_M + loc_Vb, &
                 & dgdz2_zox_M - dgdz1_zox_M + loc_Fb - w*loc_Vb /)

   loc_dim = 4
   CALL sub_matchsoln_PO4_M(loc_mat_X_4x4, loc_mat_Y_4x4, loc_vec_Z_4, loc_dim, loc_mat_C_4x4, loc_vec_D_4)

ELSE    ! 2. CASE: 3 int const. in each layer
   ! DH: zox non-bioturbated
   ! DH: this should generate 3x3 matrices as no M flux boundary condition (and then 4x4 C with zeros)
   loc_mat_X_3x3 = transpose(reshape((/ e1_zox_P, f1_zox_P, p1_zox_P, &
                          & dedz1_zox_P, dfdz1_zox_P, dpdz1_zox_P, &
                          & e1_zox_M, f1_zox_M, p1_zox_M /), shape(loc_mat_X_3x3)))

   loc_mat_Y_3x3 = transpose(reshape((/ e2_zox_P, f2_zox_P, p2_zox_P, &
                          & dedz2_zox_P, dfdz2_zox_P, dpdz2_zox_P, &
                          & e2_zox_M, f2_zox_M, p2_zox_M /), shape(loc_mat_Y_3x3)))

   loc_vec_Z_3 = (/ g2_zox_P-g1_zox_P + loc_Vb, &
                 & dgdz2_zox_P - dgdz1_zox_P + loc_Fb - w*loc_Vb, &
                 & g2_zox_M-g1_zox_M + loc_Vb /)
   loc_dim = 3

   CALL sub_matchsoln_PO4_M(loc_mat_X_3x3, loc_mat_Y_3x3, loc_vec_Z_3, loc_dim, loc_mat_C_3x3, loc_vec_D_3)

   ! integrate the 3x3 matrix in the 4x4 matrix
   DO loc_i = 1,3
      DO loc_j = 1,3
         loc_mat_C_4x4(loc_i, loc_j) = loc_mat_C_3x3(loc_i, loc_j)
      END DO
      
      loc_vec_D_4(loc_i)  = loc_vec_D_3(loc_i)
   END DO

END IF !(zox < zbio)

! Solution at SWI, top of layer 1
CALL sub_calcfg_l12_PO4_M(z=z0, reac1P=reac1_po4_ox, dum_ktempP=ksPO4/(1.0D0+KPO4_ox), &
                        & dum_QtempP=PO4s*ksPO4/(1.0D0+KPO4_ox), &
                        & dum_D1P=DPO41/(1.0D0+KPO4_ox), dum_D2P=DPO42/(1.0D0+KPO4_ox), dum_alphaP=0.0D0, &
                        & dum_mat_C=dum_mat_C1, &
                        & dum_vec_D=dum_vec_D1, dum_ltype=dum_ltype1, dum_ktempM=0.0D0, dum_QtempM=0.0D0, &
                        & dum_D1M=Dbio, &
                        & dum_D2M=0.0D0, dum_alphaM=(1.0D0/SD)*ksPO4, e_P=e1_z0_P, dedz_P=dedz1_z0_P, f_P=f1_z0_P, &
                        & dfdz_P=dfdz1_z0_P, g_P=g1_z0_P, dgdz_P=dgdz1_z0_P, p_P=p1_z0_P, dpdz_P=dpdz1_z0_P, q_P=q1_z0_P, &
                        & dqdz_P=dqdz1_z0_P, e_M=e1_z0_M, dedz_M=dedz1_z0_M, f_M=f1_z0_M, dfdz_M=dfdz1_z0_M, g_M=g1_z0_M, &
                        & dgdz_M=dgdz1_z0_M, p_M=p1_z0_M, dpdz_M=dpdz1_z0_M, q_M=q1_z0_M, dqdz_M=dqdz1_z0_M, &
                        & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
 
! transform to use coeffs from l2
! Now find 'transformed' basis functions such that in layer 1 (here for O2, PO4 is a bit more complex)
! O2 = A_2*et + B_2*ft + gt  (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
loc_EFPQ_P = (/ e1_z0_P, f1_z0_P, p1_z0_P, q1_z0_P /)
loc_dEFPQdz_P = (/ dedz1_z0_P, dfdz1_z0_P, dpdz1_z0_P, dqdz1_z0_P /)
loc_EFPQ_M = (/ e1_z0_M, f1_z0_M, p1_z0_M, q1_z0_M /)
loc_dEFPQdz_M = (/dedz1_z0_M, dfdz1_z0_M, dpdz1_z0_M, dqdz1_z0_M /)

CALL sub_xformsoln_PO4_M(loc_EFPQ_P, loc_EFPQ_M, loc_dEFPQdz_P, loc_dEFPQdz_M, &
                       & g1_z0_P, g1_z0_M, dgdz1_z0_P, dgdz1_z0_M, loc_mat_C_4x4, loc_vec_D_4, &
                       & loc_EFPQ_P_t, g_P, loc_dEFPQdz_P_t, dgdz_P, loc_EFPQ_M_t, &
                       & g_M, loc_dEFPQdz_M_t, dgdz_M)

! SD assume D2 == 0 (as q, dqdz2_zinf = 0 ) and solve for 3 unknowns
loc_mat_X_3x3 =  transpose(reshape((/ dedz2_zinf_P, dfdz2_zinf_P, dpdz2_zinf_P, &
                          & loc_EFPQ_P_t(1), loc_EFPQ_P_t(2), loc_EFPQ_P_t(3), &
                          & w*loc_EFPQ_M_t(1), w*loc_EFPQ_M_t(2), w*loc_EFPQ_M_t(3) /), shape(loc_mat_X_3x3)))
loc_vec_Z_3= (/ -dgdz2_zinf_P, dum_swiconc_PO4 - g_P, dum_swiflux_M - w*g_M  /)

! calculate the integration conctants for Layer 2
! just need it once, so actually no need for subroutine, but maybe for later
loc_dim = 3
CALL sub_solve2eqn_PO4_M(loc_mat_X_3x3, loc_vec_Z_3, rPO4_M_A2, rPO4_M_B2, rPO4_M_C2, loc_dim)

rPO4_M_D2 = 0.0
! save IC in a vector for a later calculation
loc_Layer2_IC = (/ rPO4_M_A2, rPO4_M_B2, rPO4_M_C2, rPO4_M_D2 /)

! calculate concentration at zinf
loc_conczinf = rPO4_M_A2*e2_zinf_P+rPO4_M_B2*f2_zinf_P + g2_zinf_P

! CALCULATE FINAL SWI fluxes and save the coefficients for
! DH: use A2, B2, C2, D2 as these are _xformed_ layer 1 basis functions
dum_swiflux_PO4 = por*(DPO41/(1+KPO4_ox)*(rPO4_M_A2*loc_dEFPQdz_P_t(1)+rPO4_M_B2*loc_dEFPQdz_P_t(2) + &
                & rPO4_M_C2*loc_dEFPQdz_P_t(3)+rPO4_M_D2*loc_dEFPQdz_P_t(4) + dgdz_P) - &
                & w*(dum_swiconc_PO4 - loc_conczinf))
! Does actually not exist, as it is a solid, just calculate for debugging
dum_swiflux_M = por*Dbio*(rPO4_M_A2*loc_dEFPQdz_M_t(1)+rPO4_M_B2*loc_dEFPQdz_M_t(2) + &
              & rPO4_M_C2*loc_dEFPQdz_M_t(3) + rPO4_M_D2*loc_dEFPQdz_M_t(4) + dgdz_M)

! save coeffs for layer 1 - in case I want to calculate a profile later
loc_Layer1_IC = matmul(loc_mat_C_4x4, loc_Layer2_IC) + loc_vec_D_4
!print*, loc_Layer1_IC
! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_PO4) = 0.0D0

DO j=1, OS_vertical_grid
   
   !   Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
   !   layer 1: 0 < z < zox, OM degradation (-) Sorption to sediment Fe-oxides (ktemp)
   IF (dum_z_vector(j).LE.zox) THEN 
   
   CALL sub_calcfg_l12_PO4_M(z=dum_z_vector(j), reac1P=reac1_po4_ox, dum_ktempP=ksPO4/(1.0D0+KPO4_ox), &
                           & dum_QtempP=PO4s*ksPO4/(1.0D0+KPO4_ox), &
                           & dum_D1P=DPO41/(1.0D0+KPO4_ox), dum_D2P=DPO42/(1.0D0+KPO4_ox), dum_alphaP=0.0D0, &
                           & dum_mat_C=dum_mat_C1, dum_vec_D=dum_vec_D1, dum_ltype=dum_ltype1, dum_ktempM=0.0D0, dum_QtempM=0.0D0, &
                           & dum_D1M=Dbio, dum_D2M=0.0D0, dum_alphaM=(1.0D0/SD)*ksPO4, e_P=e_zx_P, dedz_P=dedz_zx_P, &
                           & f_P=f_zx_P, &
                           & dfdz_P=dfdz_zx_P, g_P=g_zx_P, dgdz_P=dgdz_zx_P, p_P=p_zx_P, dpdz_P=dpdz_zx_P, &
                           & q_P=q_zx_P, &
                           & dqdz_P=dqdz_zx_P, e_M=e_zx_M, dedz_M=dedz_zx_M, f_M=f_zx_M, dfdz_M=dfdz_zx_M, &
                           & g_M=g_zx_M, &
                           & dgdz_M=dgdz_zx_M, p_M=p_zx_M, dpdz_M=dpdz_zx_M, q_M=q_zx_M, dqdz_M=dqdz_zx_M, &
                           & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                           
   !dum_diag_profile(j,iarr_OS_PO4) = rPO4_M_A1*e_zx_P + rPO4_M_B1*f_zx_P + rPO4_M_C1*p_zx_P + rPO4_M_D1*q_zx_P + g_zx_P
   !dum_diag_profile(j,iarr_OS_M)   = rPO4_M_A1*e_zx_M + rPO4_M_B1*f_zx_M + rPO4_M_C1*p_zx_M + rPO4_M_D1*q_zx_M + g_zx_M
   dum_diag_profile(j,iarr_OS_PO4) = loc_Layer1_IC(1)*e_zx_P + loc_Layer1_IC(2)*f_zx_P + loc_Layer1_IC(3)*p_zx_P + &
                                   & loc_Layer1_IC(4)*q_zx_P + g_zx_P
   dum_diag_profile(j,iarr_OS_M)   = loc_Layer1_IC(1)*e_zx_M + loc_Layer1_IC(2)*f_zx_M + loc_Layer1_IC(3)*p_zx_M + &
                                   & loc_Layer1_IC(4)*q_zx_M + g_zx_M
 
   !   layer 2: zox < z < zinf,
   !   OM degradation (-) authigenic P formation (ktemp) (+) P desorption due to Fe-bound P release upon Fe oxide reduction                        
   ELSEIF (dum_z_vector(j).LE.zinf) THEN

   CALL sub_calcfg_l12_PO4_M(z=dum_z_vector(j), reac1P=reac1_po4_anox, dum_ktempP=kaPO4/(1.0D0+KPO4_anox), &
                           & dum_QtempP=PO4a*kaPO4/(1.0D0+KPO4_anox), &
                           & dum_D1P=DPO41/(1.0D0+KPO4_anox), dum_D2P=DPO42/(1.0D0+KPO4_anox), &
                           & dum_alphaP=SD*kmPO4/(1.0D0+KPO4_anox), &
                           & dum_mat_C=dum_mat_C2, dum_vec_D=dum_vec_D2, dum_ltype=dum_ltype2, dum_ktempM=kmPO4, &
                           & dum_QtempM=kmPO4*Minf, &
                           & dum_D1M=Dbio, dum_D2M=0.0D0, dum_alphaM=0.0D0, e_P=e_zx_P, dedz_P=dedz_zx_P, &
                           & f_P=f_zx_P, &
                           & dfdz_P=dfdz_zx_P, g_P=g_zx_P, dgdz_P=dgdz_zx_P, p_P=p_zx_P, dpdz_P=dpdz_zx_P, &
                           & q_P=q_zx_P, &
                           & dqdz_P=dqdz_zx_P, e_M=e_zx_M, dedz_M=dedz_zx_M, f_M=f_zx_M, dfdz_M=dfdz_zx_M, &
                           & g_M=g_zx_M, &
                           & dgdz_M=dgdz_zx_M, p_M=p_zx_M, dpdz_M=dpdz_zx_M, q_M=q_zx_M, dqdz_M=dqdz_zx_M, &
                           & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                           
   !dum_diag_profile(j,iarr_OS_PO4) = rPO4_M_A2*e_zx_P + rPO4_M_B2*f_zx_P + rPO4_M_C2*p_zx_P + rPO4_M_D2*q_zx_P + g_zx_P
   !dum_diag_profile(j,iarr_OS_M)   = rPO4_M_A2*e_zx_M + rPO4_M_B2*f_zx_M + rPO4_M_C2*p_zx_M + rPO4_M_D2*q_zx_M + g_zx_M
   dum_diag_profile(j,iarr_OS_PO4) = loc_Layer2_IC(1)*e_zx_P + loc_Layer2_IC(2)*f_zx_P + loc_Layer2_IC(3)*p_zx_P + &
                                   & loc_Layer2_IC(4)*q_zx_P + g_zx_P
   dum_diag_profile(j,iarr_OS_M)   = loc_Layer2_IC(1)*e_zx_M + loc_Layer2_IC(2)*f_zx_M + loc_Layer2_IC(3)*p_zx_M + &
                                   & loc_Layer2_IC(4)*q_zx_M + g_zx_M
                        
   END IF 
      
END DO


END SUBROUTINE OMENSED_calculate_zPO4_M

!========================================================================

SUBROUTINE OMENSED_calculate_zDIC(dum_swiconc_DIC, dum_swiflux_DIC, dum_POC_conc_swi, dum_RCM_approx, &
                                & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zDIC* calculate DIC
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: dum_swiconc_DIC                ! DIC concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_DIC         ! DIC flux
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
!
REAL(8) :: loc_conczinf
REAL(8) :: reac1_dic, reac2_dic                 ! reactive terms: OM degradation
INTEGER :: ltype1 , ltype2
REAL(8) :: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
REAL(8) :: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
REAL(8) :: zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f
REAL(8) :: e2_zinf, dedz2_zinf, f2_zinf, dfdz2_zinf, g2_zinf, dgdz2_zinf
REAL(8) :: e1_zso4, dedz1_zso4, f1_zso4, dfdz1_zso4, g1_zso4, dgdz1_zso4
REAL(8) :: e2_zso4, dedz2_zso4, f2_zso4, dfdz2_zso4, g2_zso4, dgdz2_zso4
REAL(8) :: e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
REAL(8) :: e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
REAL(8) :: rDIC_A2, rDIC_B2
REAL(8) :: rDIC_A1, rDIC_B1
REAL(8) :: zso4FDIC
INTEGER :: j

reac1_dic = DICC1                 ! DIC/C until zSO4 (mol/mol)
reac2_dic = DICC2                 ! DIC/C below zSO4 (mol/mol)

print*, '---------------------- START zDIC ------------------------------- '

! Calculate DIC

! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
! layer 1: 0 < z < zso4, DIC produced my OM degradation
!    prepfg_l12(    reac1,      reac2,   ktemp, zU,  zL,   D1,    D2, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)
CALL sub_prepfg_l12(reac=reac1_dic, ktemp=0.0D0, zU=0.0D0, zL=zso4, D1=DDIC1, &
                  & D2=DDIC2, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ltype=ltype1, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! layer 2: zso4 < z < zinf, DIC production by OM degradation (Methanogenesis) -> different production rate
CALL sub_prepfg_l12(reac=reac2_dic, ktemp=0.0D0, zU=zso4, zL=zinf, D1=DDIC1, &
                  & D2=DDIC2, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ltype=ltype2, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 2 zinf
CALL sub_calcfg_l12(z=zinf, reac=reac2_dic, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, &
                  & ls_D1=DDIC1, ls_D2=DDIC2, ltype=ltype2, e=e2_zinf, dedz=dedz2_zinf, f=f2_zinf, dfdz=dfdz2_zinf, &
                  & g=g2_zinf, dgdz=dgdz2_zinf, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)


! Match at zso4, layer 1 - layer 2 (continuity and flux with AOM production)
! basis functions at bottom of layer 1
CALL sub_calcfg_l12(z=zso4, reac=reac1_dic, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, &
                  & ls_D1=DDIC1, ls_D2=DDIC2, ltype=ltype1, e=e1_zso4, dedz=dedz1_zso4, f=f1_zso4, dfdz=dfdz1_zso4, &
                  & g=g1_zso4, dgdz=dgdz1_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 2
CALL sub_calcfg_l12(z=zso4, reac=reac2_dic, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, &
                  & ls_D1=DDIC1, ls_D2=DDIC2, ltype=ltype2, e=e2_zso4, dedz=dedz2_zso4, f=f2_zso4, dfdz=dfdz2_zso4, &
                  & g=g2_zso4, dgdz=dgdz2_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! flux of DIC produced by AOM interface (Source of DIC)
IF (zso4.EQ.zinf) THEN
   zso4FDIC = 0.0D0
ELSE
   zso4FDIC = FUN_OMENSED_calcReac(zU=zso4, zL=zinf, reac=MC, dum_RCM_approx=dum_RCM_approx, &
                                 & dum_POC_conc_swi=dum_POC_conc_swi) ! MULTIPLY BY 1/POR ????!             
END IF

! match solutions at zso4 - continuous concentration and flux
CALL sub_matchsoln(e1_zso4, f1_zso4, g1_zso4, dedz1_zso4, dfdz1_zso4, dgdz1_zso4, &
                 & e2_zso4, f2_zso4, g2_zso4, dedz2_zso4, dfdz2_zso4, dgdz2_zso4, &
                 & 0.0D0, -gammaCH4*zso4FDIC/DDIC2, zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f)

! Solution at swi, top of layer 1
CALL sub_calcfg_l12(z=0.0D0, reac=reac1_dic, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, &
                  & ls_D1=DDIC1, ls_D2=DDIC2, ltype=ltype1, e=e1_00, dedz=dedz1_00, f=f1_00, dfdz=dfdz1_00, &
                  & g=g1_00, dgdz=dgdz1_00, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! transform to use coeffs from l2
CALL sub_xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
                 & zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

! Solve for ADIC, BDIC given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
!  ADIC*dedz4_zinf   +  BDIC*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
!  ADIC*e1_0     +   BDIC*f1_0     + g1_0  = swi.DIC0;
!  | dedz4_zinf dfdz4_zinf |  |ADIC|   = | -dgz4_zinf       |
!  | e1_0     f1_0         |  |BDIC|     | swi.DIC0 - g1_0 |
CALL sub_solve2eqn(dedz2_zinf, dfdz2_zinf, e1_0, f1_0, -dgdz2_zinf, dum_swiconc_DIC - g1_0, rDIC_A2, rDIC_B2)

! calculate concentration at zinf
loc_conczinf = rDIC_A2*e2_zinf+rDIC_B2*f2_zinf + g2_zinf

! flux at swi - DO include por so this is per cm^2 water column area
dum_swiflux_DIC = por*(DDIC1*(rDIC_A2*dedz1_0+rDIC_B2*dfdz1_0 + dgdz1_0) - w*(dum_swiconc_DIC - loc_conczinf))   ! NB: use A2, B2 as these are _xformed_ layer 1 basis functions

! save coeffs for layers 1
rDIC_A1 = zso4_a*rDIC_A2 + zso4_b*rDIC_B2 + zso4_e
rDIC_B1 = zso4_c*rDIC_A2 + zso4_d*rDIC_B2 + zso4_f

! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_DIC) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zso4) THEN
       
       ! layer 1: 0 < z < zso4, DIC produced my OM degradation
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac1_dic, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, &
                         & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ls_D1=DDIC1, ls_D2=DDIC2, ltype=ltype1, &
                         & e=e_zx, dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
   
       dum_diag_profile(j,iarr_OS_DIC) = rDIC_A1*e_zx+rDIC_B1*f_zx + g_zx
       
   END IF
   
   IF (dum_z_vector(j).GT.zso4) THEN
       
       ! layer 2: zso4 < z < zinf, DIC production by OM degradation (Methanogenesis) -> different production rate
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac2_dic, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                         & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DDIC1, ls_D2=DDIC2, ltype=ltype2, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_DIC) = rDIC_A2*e_zx+rDIC_B2*f_zx + g_zx
   
   END IF
      
END DO

END SUBROUTINE OMENSED_calculate_zDIC

!========================================================================

SUBROUTINE OMENSED_calculate_zALK(dum_swiconc_ALK, dum_swiflux_ALK, dum_POC_conc_swi, dum_RCM_approx, &
                                & dum_z_vector, dum_diag_profile)
!************************************************************************
!
! *OMENSED_calculate_zALK* calculate alkalinity 
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: dum_swiconc_ALK                ! ALK concentrations at SWI
REAL(8), INTENT(INOUT) :: dum_swiflux_ALK         ! ALK flux
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: dum_z_vector 
REAL(8),DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(INOUT) :: dum_diag_profile 
!
!*Local variables
!
REAL(8) :: loc_conczinf
REAL(8) :: reac11_alk, reac12_alk, reac21_alk, reac22_alk, reac3_alk, reac4_alk                 ! reactive terms: OM degradation
INTEGER :: ltype1, ltype2, ltype3, ltype4
REAL(8) :: ls_a1, ls_b1, ls_c1, ls_d1, ls_e1, ls_f1
REAL(8) :: ls_a2, ls_b2, ls_c2, ls_d2, ls_e2, ls_f2
REAL(8) :: ls_a3, ls_b3, ls_c3, ls_d3, ls_e3, ls_f3
REAL(8) :: ls_a4, ls_b4, ls_c4, ls_d4, ls_e4, ls_f4
REAL(8) :: e4_zinf, dedz4_zinf, f4_zinf, dfdz4_zinf, g4_zinf, dgdz4_zinf
REAL(8) :: e3_zso4, dedz3_zso4, f3_zso4, dfdz3_zso4, g3_zso4, dgdz3_zso4
REAL(8) :: e4_zso4, dedz4_zso4, f4_zso4, dfdz4_zso4, g4_zso4, dgdz4_zso4
REAL(8) :: zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f
REAL(8) :: e2_zno3, dedz2_zno3, f2_zno3, dfdz2_zno3, g2_zno3, dgdz2_zno3
REAL(8) :: e3_zno30, dedz3_zno30, f3_zno30, dfdz3_zno30, g3_zno30, dgdz3_zno30
REAL(8) :: e3_zno3, dedz3_zno3, f3_zno3, dfdz3_zno3, g3_zno3, dgdz3_zno3
REAL(8) :: zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f
REAL(8) :: e1_zox, dedz1_zox, f1_zox, dfdz1_zox, g1_zox, dgdz1_zox
REAL(8) :: e2_zox0, dedz2_zox0, f2_zox0, dfdz2_zox0, g2_zox0, dgdz2_zox0
REAL(8) :: e2_zox, dedz2_zox, f2_zox, dfdz2_zox, g2_zox, dgdz2_zox
REAL(8) :: zox_a, zox_b, zox_c, zox_d, zox_e, zox_f
REAL(8) :: e1_00, dedz1_00, f1_00, dfdz1_00, g1_00, dgdz1_00
REAL(8) :: e1_0, dedz1_0, f1_0, dfdz1_0, g1_0, dgdz1_0
REAL(8) :: e_zx, dedz_zx, f_zx, dfdz_zx, g_zx, dgdz_zx
REAL(8) :: rALK_A4, rALK_B4
REAL(8) :: rALK_A3, rALK_B3
REAL(8) :: rALK_A2, rALK_B2
REAL(8) :: rALK_A1, rALK_B1
REAL(8) :: zso4FALK, zoxFALK
INTEGER :: j

reac11_alk = gammaNH4*NC1/(1+KNH4)*ALKRNIT+ALKROX           ! z < zox:  Nitrification (-2) Aerobic degradation (+15/106)
reac12_alk = gammaNH4*NC2/(1+KNH4)*ALKRNIT+ALKROX           ! z < zox:  Nitrification (-2) Aerobic degradation (+15/106)
reac21_alk = ALKRDEN                                     ! zox < z < zno3: Denitrification (+93.4/106)
reac22_alk = ALKRDEN                                     ! zox < z < zno3: Denitrification (+93.4/106)
reac3_alk = ALKRSUL                                    ! zno3 < z < zso4: Sulfate reduction (+15/106)
reac4_alk = ALKRMET                                      ! zso4 < z < zinf: Methanogenesis (+14/106)

! Calculate ALK
! Preparation: for each layer, sort out solution-matching across bioturbation boundary if necessary
! layer 1: 0 < z < zox, Nitrification (-) Aerobic degradation (+)
!    prepfg_l12(reac1,      reac2,    ktemp, zU,  zL,   D1,     D2, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype)
CALL sub_prepfg_l12(reac=reac11_alk, ktemp=0.0D0, zU=0.0D0, zL=zox, D1=DALK1, &
                  & D2=DALK2, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, &
                  & ls_f=ls_f1, ltype=ltype1, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 2: zox < z < zno3, Denitrification (+)
CALL sub_prepfg_l12(reac=reac21_alk, ktemp=0.0D0, zU=zox, zL=zno3, D1=DALK1, &
                  & D2=DALK2, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, &
                  & ls_f=ls_f2, ltype=ltype2, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 3: zno3 < z < zso4, Sulfate reduction (+)
CALL sub_prepfg_l12(reac=reac3_alk, ktemp=0.0D0, zU=zno3, zL=zso4, D1=DALK1, &
                  & D2=DALK2, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, &
                  & ls_f=ls_f3, ltype=ltype3, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! layer 4: zso4 < z < zinf, Methanogenesis (+)
CALL sub_prepfg_l12(reac=reac4_alk, ktemp=0.0D0, zU=zso4, zL=zinf, D1=DALK1, &
                  & D2=DALK2, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, ls_e=ls_e4, &
                  & ls_f=ls_f4, ltype=ltype4, dum_POC_conc_swi=dum_POC_conc_swi, &
                  & dum_RCM_approx=dum_RCM_approx)

! Work up from the bottom, matching solutions at boundaries
! Basis functions at bottom of layer 4 zinf
CALL sub_calcfg_l12(z=zinf, reac=reac4_alk, ktemp=0.0D0, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, ls_e=ls_e4, ls_f=ls_f4, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype4, e=e4_zinf, dedz=dedz4_zinf, f=f4_zinf, dfdz=dfdz4_zinf, &
                  & g=g4_zinf, dgdz=dgdz4_zinf, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! Match at zso4, layer 3 - layer 4 (continuity and flux with AOM production)
! basis functions at bottom of layer 3
CALL sub_calcfg_l12(z=zso4, reac=reac3_alk, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype3, e=e3_zso4, dedz=dedz3_zso4, f=f3_zso4, dfdz=dfdz3_zso4, &
                  & g=g3_zso4, dgdz=dgdz3_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 4
CALL sub_calcfg_l12(z=zso4, reac=reac4_alk, ktemp=0.0D0, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, ls_e=ls_e4, ls_f=ls_f4, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype4, e=e4_zso4, dedz=dedz4_zso4, f=f4_zso4, dfdz=dfdz4_zso4, &
                  & g=g4_zso4, dgdz=dgdz4_zso4, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! flux of ALK produced by AOM interface (Source of ALK)
IF (zso4.EQ.zinf) THEN
   zso4FALK = 0.0D0
ELSE
   zso4FALK = ALKRAOM*gammaCH4*FUN_OMENSED_calcReac(zU=zso4, zL=zinf, reac=MC, &
                                                  & dum_RCM_approx=dum_RCM_approx, &
                                                  & dum_POC_conc_swi=dum_POC_conc_swi) ! MULTIPLY BY 1/POR ????
END IF

! match solutions at zso4 - continuous concentration and flux
CALL sub_matchsoln(e3_zso4, f3_zso4, g3_zso4, dedz3_zso4, dfdz3_zso4, dgdz3_zso4, &
                 & e4_zso4, f4_zso4, g4_zso4, dedz4_zso4, dfdz4_zso4, dgdz4_zso4, &
                 & 0.0D0, -zso4FALK/DALK2, zso4_a, zso4_b, zso4_c, zso4_d, zso4_e, zso4_f)

! Match at zno3, layer 2 - layer 3 (continuity and flux)
! basis functions at bottom of layer 2
CALL sub_calcfg_l12(z=zno3, reac=reac21_alk, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype2, e=e2_zno3, dedz=dedz2_zno3, f=f2_zno3, dfdz=dfdz2_zno3, &
                  & g=g2_zno3, dgdz=dgdz2_zno3, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! ... and top of layer 3
CALL sub_calcfg_l12(z=zno3, reac=reac3_alk, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, ls_e=ls_e3, ls_f=ls_f3, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype3, e=e3_zno30, dedz=dedz3_zno30, f=f3_zno30, dfdz=dfdz3_zno30, &
                  & g=g3_zno30, dgdz=dgdz3_zno30, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                  
! ... transformed to use coeffs from l4
CALL sub_xformsoln(e3_zno30, f3_zno30, g3_zno30, dedz3_zno30, dfdz3_zno30, dgdz3_zno30, &
                 & zso4_a , zso4_b , zso4_c , zso4_d , zso4_e ,zso4_f, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3)
! match solutions at zno3 - continuous concentration and flux
CALL sub_matchsoln(e2_zno3, f2_zno3, g2_zno3, dedz2_zno3, dfdz2_zno3, dgdz2_zno3, &
                 & e3_zno3, f3_zno3, g3_zno3, dedz3_zno3, dfdz3_zno3, dgdz3_zno3, &
                 & 0.0D0, 0.0D0, zno3_a, zno3_b, zno3_c, zno3_d, zno3_e, zno3_f)

! Match at zox, layer 1 - layer 2 (continuity, flux discontinuity from ALK source)
! flux of ALK to oxic interface (from all sources of ALK below) from NH4 and H2S
zoxFALK = ALKRNIT*gammaNH4*FUN_OMENSED_calcReac(zU=zno3, zL=zinf, reac=NC1/(1.D0+KNH4), &
                       & dum_RCM_approx=dum_RCM_approx, dum_POC_conc_swi=dum_POC_conc_swi) + &      ! was until 27.02. -1.0 before -2.0* gamma ... MULTIPLY BY 1/POR ????
        & ALKRH2S*gammaH2S*FUN_OMENSED_calcReac(zU=zno3, zL=zso4, reac=SO4C, &
                       & dum_RCM_approx=dum_RCM_approx, dum_POC_conc_swi=dum_POC_conc_swi)                 ! Dominik 25.02.2016

CALL sub_calcfg_l12(z=zox, reac=reac11_alk, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype1, e=e1_zox, dedz=dedz1_zox, f=f1_zox, dfdz=dfdz1_zox, &
                  & g=g1_zox, dgdz=dgdz1_zox, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! basis functions at top of layer 2
CALL sub_calcfg_l12(z=zox, reac=reac21_alk, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, ls_e=ls_e2, ls_f=ls_f2, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype2, e=e2_zox0, dedz=dedz2_zox0, f=f2_zox0, dfdz=dfdz2_zox0, &
                  & g=g2_zox0, dgdz=dgdz2_zox0, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

!   transform to use coeffs from l4
CALL sub_xformsoln(e2_zox0, f2_zox0, g2_zox0, dedz2_zox0, dfdz2_zox0, dgdz2_zox0, &
                 & zno3_a , zno3_b , zno3_c , zno3_d , zno3_e ,zno3_f, &
                 & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox)

! match solutions at zox - continuous concentration, flux discontinuity from ALK ox
IF (zox.LE.zbio) THEN
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, -r_zxf*zoxFALK/DALK1, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
ELSE
   CALL sub_matchsoln(e1_zox, f1_zox, g1_zox, dedz1_zox, dfdz1_zox, dgdz1_zox, &
                    & e2_zox, f2_zox, g2_zox, dedz2_zox, dfdz2_zox, dgdz2_zox, &
                    & 0.0D0, -r_zxf*zoxFALK/DALK2, zox_a, zox_b, zox_c, zox_d, zox_e, zox_f)
END IF

! Solution at swi, top of layer 1
CALL sub_calcfg_l12(z=0.0D0, reac=reac11_alk, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, &
                  & ls_D1=DALK1, ls_D2=DALK2, ltype=ltype1, e=e1_00, dedz=dedz1_00, f=f1_00, dfdz=dfdz1_00, &
                  & g=g1_00, dgdz=dgdz1_00, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

! transform to use coeffs from l4
CALL sub_xformsoln(e1_00, f1_00, g1_00, dedz1_00, dfdz1_00, dgdz1_00, &
                 & zox_a , zox_b , zox_c , zox_d , zox_e ,zox_f, &
                 & e1_0, f1_0, g1_0, dedz1_0,  dfdz1_0, dgdz1_0)

! Solve for AALK, BALK given boundary conditions (expressed in terms of transformed basis fns, layer 4 A, B)
!  AALK*dedz4_zinf   +  BALK*dfz4_zinf  + dgz4_zinf = 0;          % zero flux at zinf
!  AALK*e1_0     +   BALK*f1_0     + g1_0  = swi.ALK0;
!  | dedz4_zinf dfdz4_zinf |  |AALK|   = | -dgz4_zinf       |
!  | e1_0     f1_0         |  |BALK|     | swi.ALK0 - g1_0 |

CALL sub_solve2eqn(dedz4_zinf, dfdz4_zinf, e1_0, f1_0, -dgdz4_zinf, dum_swiconc_ALK - g1_0, rALK_A4, rALK_B4)

! calculate concentration at zinf
loc_conczinf = rALK_A4*e4_zinf+rALK_B4*f4_zinf + g4_zinf

! flux at swi - DO include por so this is per cm^2 water column area
dum_swiflux_ALK = por*(DALK1*(rALK_A4*dedz1_0+rALK_B4*dfdz1_0 + dgdz1_0) - w*(dum_swiconc_ALK - loc_conczinf))   ! NB: use A4, B4 as these are _xformed_ layer 1 basis functions

! save coeffs for layers 3, 2 and 1
rALK_A3 = zso4_a*rALK_A4 + zso4_b*rALK_B4 + zso4_e
rALK_B3 = zso4_c*rALK_A4 + zso4_d*rALK_B4 + zso4_f
rALK_A2 = zno3_a*rALK_A4 + zno3_b*rALK_B4 + zno3_e
rALK_B2 = zno3_c*rALK_A4 + zno3_d*rALK_B4 + zno3_f
rALK_A1 = zox_a*rALK_A4 + zox_b*rALK_B4 + zox_e
rALK_B1 = zox_c*rALK_A4 + zox_d*rALK_B4 + zox_f

! Fill in diagenetic profile array for later plotting
dum_diag_profile(:,iarr_OS_AT) = 0.0D0

DO j=1, OS_vertical_grid
   
   IF (dum_z_vector(j).LE.zox) THEN
       
       ! layer 1: 0 < z < zox, Nitrification (-) Aerobic degradation (+)
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac11_alk, ktemp=0.0D0, ls_a=ls_a1, ls_b=ls_b1, &
                         & ls_c=ls_c1, ls_d=ls_d1, ls_e=ls_e1, ls_f=ls_f1, ls_D1=DALK1, ls_D2=DALK2, ltype=ltype1, &
                         & e=e_zx, dedz=dedz_zx, f=f_zx, dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
   
       dum_diag_profile(j,iarr_OS_AT) = rALK_A1*e_zx+rALK_B1*f_zx + g_zx
       
   ELSEIF (dum_z_vector(j) .LE. zno3) THEN !dum_z_vector(j).GT.zox .AND. 
       
       ! layer 2: zox < z < zno3, Denitrification (+)
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac21_alk, ktemp=0.0D0, ls_a=ls_a2, ls_b=ls_b2, ls_c=ls_c2, ls_d=ls_d2, &
                         & ls_e=ls_e2, ls_f=ls_f2, ls_D1=DALK1, ls_D2=DALK2, ltype=ltype2, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_AT) = rALK_A2*e_zx+rALK_B2*f_zx + g_zx
   
   ELSEIF (dum_z_vector(j) .LE. zso4) THEN 
       
       ! layer 3: zno3 < z < zso4, Sulfate reduction (+)
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac3_alk, ktemp=0.0D0, ls_a=ls_a3, ls_b=ls_b3, ls_c=ls_c3, ls_d=ls_d3, &
                         & ls_e=ls_e3, ls_f=ls_f3, ls_D1=DALK1, ls_D2=DALK2, ltype=ltype3, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_AT) = rALK_A3*e_zx+rALK_B3*f_zx + g_zx

   ELSE        
       ! layer 4: zso4 < z < zinf, Methanogenesis (+)
       CALL sub_calcfg_l12(z=dum_z_vector(j), reac=reac4_alk, ktemp=0.0D0, ls_a=ls_a4, ls_b=ls_b4, ls_c=ls_c4, ls_d=ls_d4, &
                         & ls_e=ls_e4, ls_f=ls_f4, ls_D1=DALK1, ls_D2=DALK2, ltype=ltype4, e=e_zx, dedz=dedz_zx, f=f_zx, &
                         & dfdz=dfdz_zx, g=g_zx, dgdz=dgdz_zx, &
                         & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)

       dum_diag_profile(j,iarr_OS_AT) = rALK_A4*e_zx+rALK_B4*f_zx + g_zx
   
   END IF
   
END DO

END SUBROUTINE OMENSED_calculate_zALK

!========================================================================

FUNCTION FUN_OMENSED_calcReac(zU, zL, reac, dum_RCM_approx, dum_POC_conc_swi)
!************************************************************************
!
! *FUN_OMENSED_calcReac* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: zU, zL, reac
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAl(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcReac

! Integral of reacted organic matter from zU to zL,
! multiplied by stoichiometric factors reac1, reac2 (for the two OC phases)

! Vector-friendly way of handling 3 cases:
! 1) wholly within bioturbated layer:    FUN_OMENSED_calcReac_l1(zU,zL)     + (0 =) FUN_OMENSED_calcReac_l2(zbio, zbio)
! 2) wholly within non-bio     layer:  (0=) FUN_OMENSED_calcReac_l1(zbio, zbio) +   FUN_OMENSED_calcReac_l2(zU, zL)
! 3) crossing zbio                       calcRead_l1(zU,zbio)   +       FUN_OMENSED_calcReac_l2(zbio, zL)

FUN_OMENSED_calcReac = FUN_OMENSED_calcReac_l1(zU=min(zU,zbio), zL=min(zL,zbio), reac=reac, &
                     & dum_RCM_approx=dum_RCM_approx, dum_POC_conc_swi=dum_POC_conc_swi) + &
                     & FUN_OMENSED_calcReac_l2(zU=max(zU,zbio), zL=max(zL,zbio), reac=reac, &
                     & dum_RCM_approx=dum_RCM_approx)
! TODO confirm (1-por)*  has been added (to k1 & k2 ?)

END FUNCTION FUN_OMENSED_calcReac

!========================================================================

FUNCTION FUN_OMENSED_calcReac_l1(zU, zL, reac, dum_RCM_approx, dum_POC_conc_swi)
!************************************************************************
!
! *FUN_OMENSED_calcReac_l1* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_calcReac
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: zU, zL, reac
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
REAl(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcReac_l1
REAL(8), DIMENSION(OS_RCM_fracs) :: reacf, Aux_calcReac_l1
INTEGER :: i

DO i=1, OS_RCM_fracs
   
   reacf(i) = dum_RCM_approx(2,i)*reac
   
   Aux_calcReac_l1(i) = -reacf(i)*(A11(i)*(DEXP(aa1(i)*zU)*bb1(i) - DEXP(bb1(i)*zU)*aa1(i) - DEXP(aa1(i)*zL)*bb1(i) + &
                      & DEXP(bb1(i)*zL)*aa1(i)) + dum_POC_conc_swi(i)*DEXP(bb1(i)*zU)*aa1(i) - &
                      & dum_POC_conc_swi(i)*DEXP(bb1(i)*zL)*aa1(i))/(aa1(i)*bb1(i) + const_real_nullsmall) 
   
END DO


FUN_OMENSED_calcReac_l1 = SUM(Aux_calcReac_l1)

END FUNCTION FUN_OMENSED_calcReac_l1

!========================================================================

FUNCTION FUN_OMENSED_calcReac_l2(zU, zL, reac, dum_RCM_approx)
!************************************************************************
!
! *FUN_OMENSED_calcReac_l2* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_calcReac
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: zU, zL, reac
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcReac_l2
REAL(8), DIMENSION(OS_RCM_fracs) :: reacf, Aux_calcReac_l2

INTEGER :: i

DO i=1, OS_RCM_fracs
   
   reacf(i) = dum_RCM_approx(2,i)*reac
   
   Aux_calcReac_l2(i) = -reacf(i)*A22(i)*(DEXP(aa2(i)*zU) - DEXP(aa2(i)*zL))/(aa2(i) + const_real_nullsmall) 
      
END DO

FUN_OMENSED_calcReac_l2 = SUM(Aux_calcReac_l2)
  
END FUNCTION FUN_OMENSED_calcReac_l2

!========================================================================

FUNCTION FUN_OMENSED_calcOM(zU, zL, reac, dum_POC_conc_swi)
!************************************************************************
!
! *FUN_OMENSED_calcOM* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: zU, zL, reac
REAl(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcOM
REAL(8) :: loc_FUN_calcOM_l1, loc_FUN_calcOM_l2

! Integral of organic matter from zU to zL,
! multiplied by stoichiometric factors reac1, reac2 (for the two OC phases)

! Vector-friendly way of handling 3 cases:
! 1) wholly within bioturbated layer:    FUN_OMENSED_calcReac_l1(zU,zL)     + (0 =) FUN_OMENSED_calcReac_l2(zbio, zbio)
! 2) wholly within non-bio     layer:  (0=) FUN_OMENSED_calcReac_l1(zbio, zbio) +   FUN_OMENSED_calcReac_l2(zU, zL)
! 3) crossing zbio                       calcRead_l1(zU,zbio)   +       FUN_OMENSED_calcReac_l2(zbio, zL)

loc_FUN_calcOM_l1 = FUN_OMENSED_calcOM_l1(zU=min(zU,zbio), zL=min(zL,zbio), reac=reac, &
                                        & dum_POC_conc_swi=dum_POC_conc_swi)
loc_FUN_calcOM_l2 = FUN_OMENSED_calcOM_l2(zU=max(zU,zbio), zL=max(zL,zbio), reac=reac)

IF (ISNAN(loc_FUN_calcOM_l1).AND.ISNAN(loc_FUN_calcOM_l2)) THEN
   
   FUN_OMENSED_calcOM = 0.0D0

ELSEIF (ISNAN(loc_FUN_calcOM_l2)) THEN
   FUN_OMENSED_calcOM = loc_FUN_calcOM_l1
ELSE
   FUN_OMENSED_calcOM = loc_FUN_calcOM_l1 + loc_FUN_calcOM_l2
END IF

!TODO confirm (1-por)*  has been added (to k1 & k2 ?)

END FUNCTION FUN_OMENSED_calcOM

!========================================================================

FUNCTION FUN_OMENSED_calcOM_l1(zU, zL, reac, dum_POC_conc_swi)
!************************************************************************
!
! *FUN_OMENSED_calcOM_l1* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_calcOM
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: zU, zL, reac
REAl(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcOM_l1, reacf
REAL(8), DIMENSION(OS_RCM_fracs) :: Aux_calcOM_l1
INTEGER :: i

reacf = reac

DO i=1, OS_RCM_fracs
    
    Aux_calcOM_l1(i) = -reacf*(A11(i)*(DEXP(aa1(i)*zU)*bb1(i) - DEXP(bb1(i)*zU)*aa1(i) - &
                     & DEXP(aa1(i)*zL)*bb1(i) + DEXP(bb1(i)*zL)*aa1(i)) + &
                     & dum_POC1_conc_swi*DEXP(bb1(i)*zU)*aa1(i) - dum_POC1_conc_swi*DEXP(bb1(i)*zL)*&
                     & aa1(i))/(aa1(i)*bb1(i) + const_real_nullsmall) 
                     
END DO

FUN_OMENSED_calcOM_l1 = SUM(Aux_calcOM_l1)

END FUNCTION FUN_OMENSED_calcOM_l1

!========================================================================

FUNCTION FUN_OMENSED_calcOM_l2(zU, zL, reac)
!************************************************************************
!
! *FUN_OMENSED_calcOM_l2* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - FUN_OMENSED_calcOM
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: zU, zL, reac
!
!*Local variables
!
REAL(8) :: FUN_OMENSED_calcOM_l2, reacf
REAL(8), DIMENSION(OS_RCM_fracs) :: Aux_calcOM_l2
INTEGER :: i

reacf = reac

DO i=1, OS_RCM_fracs
    
    Aux_calcOM_l2(i) = -reacf*A22(i)*(DEXP(aa2(i)*zU) - DEXP(aa2(i)*zL))/ &
                    & (aa2(i) + const_real_nullsmall)
                     
END DO

FUN_OMENSED_calcOM_l2 = SUM(Aux_calcOM_l2)

END FUNCTION FUN_OMENSED_calcOM_l2

!========================================================================

END MODULE OMENSED_module

