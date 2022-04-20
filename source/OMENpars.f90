MODULE OMENpars
!************************************************************************
!
! *CGEMpars* CGEM estuarine model parameters
!
! Author - Sebastiaan van de Velde
!
! Last update - 20 Jan 2021  @(CGEM)CGEMpars.f90  V0.1
!
! Description -
!
!************************************************************************
!

IMPLICIT NONE

!---main parameters (from CGEM)

INTEGER :: MAXT, DELTI, nt
CHARACTER*256, PARAMETER :: logdir = "CGEM_log/"
REAL(8) :: water_temp

!---OMENSED global parameters

REAL(8) :: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f

INTEGER :: MaxOMENSWIArids, MaxOMENbcArids, OS_RCM_fracs
INTEGER :: OS_vertical_grid
REAL(8) :: RCM_a, RCM_nu

REAL(8) :: const_real_nullsmall, const_R_SI, OS_tolcte   
REAL(8) :: conv_POC_cm3_mol, fun_calc_sed_vol
REAL(8) :: par_sed_huelse2017_k1, par_sed_huelse2017_k2, &
         & par_sed_huelse2017_k2_anoxic, par_sed_huelse2017_k2_order
CHARACTER :: par_sed_huelse2017_kscheme
LOGICAL :: par_sed_huelse2017_P_cycle, par_sed_huelse2017_redox, &
         & par_sed_huelse2017_remove_impl_sulalk, par_sed_huelse2017_sim_P_loss, &
         & par_sed_huelse2017_sim_P_loss_pres_fracC, par_sed_huelse2017_sim_P_regeneration

LOGICAL :: OS_remineralize_all      ! remineralize everything at SWI, i.e. reflective boundary
!
! Name        Type    Purpose                                                 Unit
!------------------------------------------------------------------------------
!*TOL*        REAL    convergence criterium
!
!*nt*         INTEGER timestep                                                [s]
!*MAXT*       INTEGER max time + add hydrodynamic & transport warmup days     [s]
!*WARMUP*     INTEGER warm up period; switch on biogeochemistry if t>WARMUP   [s]
!*DELTI*      INTEGER delta t                                                 [s]
!*TS*         INTEGER save every TS timestep (1= every 150 seconds; 4= every 10 minutes; 8=ever 20 minutes; 96=every 4 hours)
!*MAXITS*     INTEGER max number of iteration steps

!---OMENSED sediment characteristics parameters

REAL(8) :: rho_sed, wdepth, w, z0, zox, zbio, zinf, Dbio, Dunbio, por, tort, &
         & irrigationFactor, dispFactor
REAL(8) ::  X_C, Y_N, Z_P, SD, OC, NC1, NC2, PC1, PC2, SO4C, O2H2S, DICC1, DICC2, &
         &  MC, gammaNH4, gammaH2S, gammaCH4, satSO4, NO3CR                                  
REAL(8) ::  ALKROX, ALKRNIT, ALKRDEN, ALKRSUL, ALKRH2S, ALKRMET, ALKRAOM, zoxgf   
REAL(8) ::  dum_POC1_conc_swi, dum_POC2_conc_swi, dum_POC3_conc_swi , dum_POC1_conc_swi_nonbio, &
         &  dum_POC2_conc_swi_nonbio, dum_POC_total_flux_zinf, dum_swiconc_O2, dum_swiconc_SO4, &
         &  dum_swiconc_H2S, dum_swiconc_NO3, dum_swiconc_NH4, dum_swiconc_PO4, dum_swiflux_M, &
         &  dum_swiconc_DIC, dum_swiconc_ALK              
REAL(8) :: DC1, DC2, k1, k2, qdispO2, adispO2, DO21, DO22, r_zxf                         
REAL(8) :: qdispNO3, adispNO3, DN1, DN2, zno3, KNH4, &
         & qdispSO4, adispSO4, DSO41, DSO42, zso4, &
         & qdispNH4, adispNH4, DNH41, DNH42, &
         & qdispH2S, adispH2S, DH2S1, DH2S2, &
         & qdispPO4, adispPO4, DPO41, DPO42, KPO4_ox, KPO4_anox, ksPO4, kmPO4, kaPO4, PO4s, PO4a, Minf, &
         & qdispDIC, adispDIC, DDIC1, DDIC2, &
         & qdispALK, adispALK, DALK1, DALK2 
REAL(8) :: OS_POC_remin_K1, OS_POC_remin_K2   
REAL(8) :: OS_POC_remin_Ea1, OS_POC_remin_Ea2  
REAL(8) :: OS_BW_O2_anoxia           

LOGICAL :: print_results = .FALSE.
    
!
! Name        Type    Purpose                                                 Unit
!------------------------------------------------------------------------------
!*rho_sed*    REAL    sediment density                                       [g/cm3]
!*wdepth*     REAL    water depth                                            [m]
!*w*          REAL    burial velocity                                        [cm/yr]
!*z0*         REAL    surface                                                [cm]
!*zox*        REAL    oxygen penetration depth                               [cm]
!*zbio*       REAL    bioturbation depth                                     [cm]
!*zinf*       REAL    Infinite depth                                         [cm]
!*Dbio*       REAL    bioturbation coefficient                               [cm2/yr]
!*Dunbio*     REAL    2nd diffusion coefficient in case used for 2nd Layer   [cm2/yr]
!*por*        REAL    porosity                                               [-]
!*tort*       REAL    tortuosity                                             [-]
!*irrigationFactor* REAL irrigation factor                                   [-]
!*dispFactor* REAL    dispersion factor                                      [-]
!*X_C*        REAL    Redfield ratio (POC:PON:POP)                           [-]
!*Y_N*        REAL    Redfield ratio (POC:PON:POP)                           [-]
!*Z_P*        REAL    Redfield ratio (POC:PON:POP)                           [-]
!*SD*         REAL    volume factor solid->dissolved phase                   [-]
!*OC*         REAL    O2/C                                                   [mol/mol]
!*NC1*        REAL    N/C first TOC fraction                                 [mol/mol]
!*NC2*        REAL    N/C second TOC fraction                                [mol/mol]
!*PC1*        REAL    P/C first TOC fraction                                 [mol/mol]
!*PC2*        REAL    P/C second TOC fraction                                [mol/mol]
!*SO4C*       REAL    SO4/C                                                  [mol/mol]
!*O2H2S*      REAL    O2/H2S O2 needed to oxidize H2S                        [mol/mol] 
!*DICC1*      REAL    DIC/C until zSO4                                       [mol/mol]
!*DICC2*      REAL    DIC/C below zSO4                                       [mol/mol]
!*MC*         REAL    CH4/C                                                  [mol/mol]
!*gammaNH4*   REAL    fraction of NH4 that is oxidised in oxic layer         [-]
!*gammaH2S*   REAL    fraction of H2S that is oxidised in oxic layer         [-]
!*gammaCH4*   REAL    fraction of CH4 that is oxidised at SO4                [-]
!*satSO4*     REAL    SO4 saturation                                         [mol cm-3]
!*NO3CR*      REAL    NO3 consumed by Denitrification                        [?]
!*ALKROX*     REAL    Alk prod/con during Aerobic degradation                [mol/mol]
!*ALKRNIT*    REAL    Alk prod/con during Nitrification                      [mol/mol]
!*ALKRDEN*    REAL    Alk prod/con during Denitrification                    [mol/mol]
!*ALKRSUL*    REAL    Alk prod/con during Sulfate reduction                  [mol/mol]
!*ALKRH2S*    REAL    Alk prod/con during H2S oxidation (CHECK THIS VALUE!!!) [mol/mol]
!*ALKRMET*    REAL    Alk prod/con during Methanogenesis                     [mol/mol]
!*ALKRAOM*    REAL    Alk prod/con during AOM                                [mol/mol]
!*zoxgf*      REAL    rolloff NH4, H2S oxidation for small zox depth         [cm]
!*dum_POC1_conc_swi* REAL conversion TOC flux at SWI                         (mol/(cm2 yr)) -> (mol/cm3 bulk phase)
!*dum_POC2_conc_swi* REAL conversion TOC flux at SWI                         (mol/(cm2 yr)) -> (mol/cm3 bulk phase)
!*dum_POC3_conc_swi* REAL conversion TOC flux at SWI                         (mol/(cm2 yr)) -> (mol/cm3 bulk phase)
!*dum_POC1_conc_swi_nonbio* REAL conversion TOC flux at SWI                  (mol/(cm2 yr)) -> (mol/cm3 bulk phase) not taking into account biodiffusion
!*dum_POC2_conc_swi_nonbio* REAL conversion TOC flux at SWI                  (mol/(cm2 yr)) -> (mol/cm3 bulk phase) not taking into account biodiffusion 
!*dum_POC_total_flux_zinf* REAL conversion total TOC flux at SWI             (mol/(cm2 yr)) -> (mol/cm3 bulk phase)
!*dum_swiconc_O2*  REAL  O2 concentration at SWI                             [mol/cm3]
!*dum_swiconc_SO4* REAL  SO4 concentration at SWI                            [mol/cm3]
!*dum_swiconc_H2S* REAL  H2S concentration at SWI                            [mol/cm3]
!*dum_swiconc_NO3* REAL  NO3 concentration at SWI                            [mol/cm3]
!*dum_swiconc_NH4* REAL  NH4 concentration at SWI                            [mol/cm3]
!*dum_swiconc_PO4* REAL  PO4 concentration at SWI                            [mol/cm3]
!*dum_swiflux_M*   REAL  conversion flux of M to the sediment                (mol/(cm2*yr)) -> is converted into concentration (mol/cm3)
!*dum_swiconc_DIC* REAL  DIC concentration at SWI                            [mol/cm3]
!*dum_swiconc_ALK* REAL  ALK concentration at SWI                            [mol/cm3] 
!*DC1*        REAL       TOC diffusion coefficient                           [cm2/yr]
!*DC2*        REAL       TOC diffusion coefficient                           [cm2/yr]
!*k1*         REAL       TOC degradation rate constnat                       [1/yr]
!*k2*         REAL       TOC degradation rate constnat                       [1/yr]
!*qdispO2*    REAL       O2 diffusion coefficient                            [cm2/yr]
!*adispO2*    REAL       O2 linear coefficient for temperature dependence    [cm2/yr/deg C]
!*DO21*       REAL       O2 diffusion coefficient in bioturbated layer       [cm2/yr]
!*DO22*       REAL       O2 diffusion coefficient in non-bioturbated layer   [cm2/yr]
!*r_zxf*      REAL       roll off oxidation at low zox
!*qdispNO3*   REAL       NO3 diffusion coefficient in water                  [cm2/yr]
!*adispNO3*   REAL       NO3 linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DN1*        REAL       NO3 diffusion coefficient in bioturbated layer      [cm2/yr]
!*DN2*        REAL       NO3 diffusion coefficient in non-bioturbated layer  [cm2/yr]
!*zno3*       REAL       NO3 penetration depth                               [cm]               
!*KNH4*       REAL       Adsorption coefficient (same in ocix and anoxic layer) [-]
!*qdispSO4*   REAL       SO4 diffusion coefficient in water                  [cm2/yr]
!*adispSO4*   REAL       SO4 linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DSO41*      REAL       SO4 diffusion coefficient in bioturbated layer      [cm2/yr] 
!*DSO42*      REAL       SO4 diffusion coefficient in non-bioturbated layer  [cm2/yr]
!*zso4*       REAL       SO4 penetration depth                               [cm] 
!*qdispNH4*   REAL       NH4 diffusion coefficient in water                  [cm2/yr] 
!*adispNH4*   REAL       NH4 linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DNH41*      REAL       NH4 diffusion coefficient in bioturbated layer      [cm2/yr] 
!*DNH42*      REAL       NH4 diffusion coefficient in non-bioturbated layer  [cm2/yr] 
!*qdispH2S*   REAL       H2S diffusion coefficient in water                  [cm2/yr] 
!*adispH2S*   REAL       H2S linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DH2S1*      REAL       H2S diffusion coefficient in bioturbated layer      [cm2/yr] 
!*DH2S2*      REAL       H2S diffusion coefficient in non-bioturbated layer  [cm2/yr] 
!*qdispPO4*   REAL       PO4 diffusion coefficient in water                  [cm2/yr] 
!*adispPO4*   REAL       PO4 linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DPO41*      REAL       PO4 diffusion coefficient in bioturbated layer      [cm2/yr] 
!*DPO42*      REAL       PO4 diffusion coefficient in non-bioturbated layer  [cm2/yr] 
!*KPO4_ox*    REAL       Adsorption coefficient in oxic layer                [-]
!*KPO4_anox*  REAL       Adsorption coefficient in anoxic layer              [-]
!*ksPO4*      REAL       Rate constant for kinetic P sorption                [1/yr]
!*kmPO4*      REAL       Rate constant for Fe-bound P release upon Fe oxide reduction [1/yr]
!*kaPO4*      REAL       Rate constant for authigenic P formation            [1/yr]
!*PO4s*       REAL       Equilibrium concentration for P sorption            [mol/cm3]
!*PO4a*       REAL       Equilibrium concentration for authigenic P formation[mol/cm3]
!*Minf*       REAL       asymptotic concentration for Fe-bound P             [mol/cm3]
!*qdispDIC*   REAL       DIC diffusion coefficient in water                  [cm2/yr] 
!*adispDIC*   REAL       DIC linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DDIC1*      REAL       DIC diffusion coefficient in bioturbated layer      [cm2/yr] 
!*DDIC2*      REAL       DIC diffusion coefficient in non-bioturbated layer  [cm2/yr] 
!*qdispALK*   REAL       ALK diffusion coefficient in water                  [cm2/yr] 
!*adispALK*   REAL       ALK linear coefficient for temperature dependence   [cm2/yr/deg C]
!*DALK1*      REAL       ALK diffusion coefficient in bioturbated layer      [cm2/yr] 
!*DALK2*      REAL       ALK diffusion coefficient in non-bioturbated layer  [cm2/yr] 
!*por*        REAL       porosity (porewater_vol./(solid_sed_vol.+porewater_vol.)) [-]
!*OS_POC_remin_K1*  REAL  rate constants for temperature dependent degradation
!*OS_POC_remin_K2*  REAL  rate constants for temperature dependent degradation
!*OS_POC_remin_Ea1* REAL activation energies for temperature dependent degradation
!*OS_POC_remin_Ea2* REAL activation energies for temperature dependent degradation
!*OS_BW_O2_anoxia*  REAL  BW [O2} threshold for zbio switch to 0.01cm         [mol cm-3]
!*print_results* LOGICAL Print output?                                        -
    
SAVE

END MODULE OMENpars
