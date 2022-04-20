MODULE default_OMENSED
!************************************************************************
!
! *default_OMENSED* Default settings for the OMENSED sediment model
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Last update - 20 Jan 2021  @(OMENSED)default_OMENSED.f90  V0.1
!
! Description -
!
! Routines - default_OMENSED_params
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE default_OMENSED_params
!************************************************************************
!
! *default_OMENSED_params* Default settings for OMENSED sediment model
!
! Author - Sebastiaan van de Velde
!
! Version - @(OMENSED)default_OMENSED.f90  V0.1
!
! Description -
!
! Calling program - initialise model
!
!************************************************************************
!
USE OMENpars

IMPLICIT NONE

!
!1. Model parameters
!-------------------
!
!---from CGEM
MAXT=4
DELTI=1

MaxOMENSWIArids=11
MaxOMENbcArids=9
OS_vertical_grid=100
water_temp=273.15D0 + 20.D0 ! K

!---RCM parameters
OS_RCM_fracs=100
RCM_a=1.D0
RCM_nu=0.125D0

!---OMENSED global parameters

OS_tolcte = 1D-18
const_real_nullsmall = 1D-24
const_R_SI = 0.   
conv_POC_cm3_mol = 1.
fun_calc_sed_vol = 1.

par_sed_huelse2017_P_cycle = .TRUE.
par_sed_huelse2017_redox = .FALSE.
par_sed_huelse2017_remove_impl_sulalk = .FALSE.
par_sed_huelse2017_sim_P_loss = .FALSE.
par_sed_huelse2017_sim_P_loss_pres_fracC = .FALSE.
par_sed_huelse2017_sim_P_regeneration = .FALSE.
OS_remineralize_all  = .FALSE.

!---OMENSED global parameters

rho_sed = 2.6D0                                      ! sediment density (g/cm3)
z0 = 0.0D0                                           ! top of the sediments
zox = 0.0D0
zno3 = 0.0D0
zso4 = 0.0D0
zinf = 100.0D0                                        ! Inifinity - bottom of the sediments (cm)
zbio = 1.0D0                                          ! bioturbation depth (cm)
Dbio = 0.02D0                                         !bioturbation coefficient (cm2/yr) - after Middelburg at al. 1997
Dunbio = 0.01D0                                       ! "diffusion coefficient" for unbioturbated layer - in case we use it
por= 0.85D0
tort = 3.0D0                                          ! tortuosity (-)
irrigationFactor = 1.0D0
gammaNH4 = 1.0D0 !0.9D0                                      ! fraction of NH4 that is oxidised in oxic layer
gammaH2S = 1.0D0 !0.95D0                                     ! fraction of H2S that is oxidised in oxic layer
gammaCH4 = 1.0D0     !0.99                            ! fraction of CH4 that is oxidised at SO4
satSO4 = 0.0D0                                        ! SO4 saturation
KNH4 = 1.3D0                                          ! Adsorption coefficient (same in oxic and anoxic layer) (-)
zoxgf = 0.1D0                                         ! cm, rolloff NH4, H2S oxidation for small zox depth
r_zxf=0.0D0
dispFactor=por**(tort-1.0D0)*irrigationFactor         ! dispersion factor (-) - Ausbreitung - type of mixing that accompanies 
                                                      ! hydrodynamic flows -> ~builds out paths of flow
SD=(1.0D0-por)/por                                    ! volume factor solid->dissolved phase
OC= SD*1.0D0!(138.0D0/106.0D0)!1.0*SD                       ! O2/C (mol/mol)
NC1=16.0/106.0*SD  ! 0.1509*SD                 ! N/C first TOC fraction (mol/mol)
NC2=16.0/106.0*SD  !0.13333*SD                 ! N/C second TOC fraction (mol/mol)
PC1=SD*1.0D0/106.0D0 !0.0094*SD                       ! P/C first TOC fraction (mol/mol)
PC2=SD*1.0D0/106.0D0 !0.0094*SD                       ! P/C second TOC fraction (mol/mol)
SO4C=SD*0.5D0!(138.0D0/212.0D0)!0.5*SD                      ! SO4/C (mol/mol)
O2H2S=2.0D0                                           ! 2 mole O2 oxidize 1 mole H2S
DICC1=1.0D0*SD                                        ! DIC/C until zSO4 (mol/mol)
DICC2=0.5D0*SD                                        ! DIC/C below zSO4 (mol/mol)
MC=0.5D0*SD                                           ! CH4/C (mol/mol)
NO3CR=SD*4.0D0/5.0D0!(94.4D0/106.0D0)*SD                             ! NO3 consumed by Denitrification
X_C=106.0D0                                           ! Carbon Redfield stoichiometry
Y_N=16.0D0                                            ! Nitrogen Redfield stoichiometry
Z_P=1.0D0                                             ! Phosphorous Redfield stoichiometry
ALKROX= -((Y_N)/X_C)*SD                               ! Aerobic degradation -16/106*SD
ALKRNIT=0.0D0                                         ! Nitrification explicit -2.0
ALKRDEN=0.0D0                                         ! Denitrification explicit: (4*X_C+3*Y_N-10*Z_P)/(5*X_C)*SD
ALKRSUL= ((X_C+Y_N)/X_C)*SD                           ! ((X_C+Y_N)/X_C)*SD = +122/106*SD!,  Sulfate reduction 
                                                      ! (N explicit: ((X_C+Y_N-2*Z_P)/X_C)*SD = +120/106*SD)
ALKRH2S= -2.0D0                                       ! H2S oxidation
ALKRMET= -((Y_N)/X_C)*SD                              ! Methanogenesis explicitly: ((Y_N-2*Z_P)/X_C)*SD
ALKRAOM= 2.0D0                                        ! AOM
DC1 = Dbio
DC2 = Dunbio

qdispO2=348.62172D0
adispO2=14.08608D0
qdispNO3=308.42208D0
adispNO3=12.2640D0
qdispNH4=309.0528D0
adispNH4=12.2640D0
qdispSO4=157.68D0                                     ! SO4 diffusion coefficient in water at 0 degree C  (cm2/yr)
adispSO4=7.884D0                                      ! SO4 linear coefficient for temperature dependence (cm2/yr/oC)
qdispH2S=307.476D0
adispH2S=9.636D0
qdispPO4=112.90764D0
adispPO4=5.586252D0
qdispDIC=151.69D0                                     ! DIC diffusion coefficient in water (cm2/yr)
adispDIC=7.93D0                                       ! DIC linear coefficient for temperature dependence (cm2/yr/oC)
qdispALK=151.69D0                                     ! ALK diffusion coefficient in water (cm2/yr)
adispALK=7.93D0                                       ! ALK linear coefficient for temperature dependence (cm2/yr/oC)

OS_BW_O2_anoxia = 5.0D-9
OS_POC_remin_K1 = 9.0D11
OS_POC_remin_K2 = 1.0D14
OS_POC_remin_Ea1 = 55000.0D0
OS_POC_remin_Ea2 = 80000.0D0

!loc_TempC = dum_TempK - 273.15D0
KPO4_ox = 200.0D0                                      ! Adsorption coefficient in oxic layer (-)
KPO4_anox = 1.3D0                                      ! Adsorption coefficient in anoxic layer (-)
dum_swiflux_M = 365.0D0*0.2D-10                        ! Flux input 365*0.2e-10 flux of M to the sediment (mol/(cm2*yr)) 
kmPO4 = 0.19D0                                         ! Rate constant for Fe-bound P release upon Fe oxide reduction 
                                                       ! (from Slomp et al. (1996)
kaPO4 = 0.37D0                                         ! Rate constant for authigenic P formation (1/yr)
PO4a = 3.7D-9     ! 0.0                                ! Equilibrium concentration for authigenic P formation (mol/cm3)
!IF(WaterDepth  .LE. 2000.0D0) THEN                           ! sediment margins
   ksPO4 = 36.5D0                                      ! Rate constant for kinetic P sorption (1/yr)
   PO4s = 2.0D-9                                       ! Equilibrium concentration for P sorption (mol/cm3)
   Minf = 5.2D-15                                      ! asymptotic concentration for Fe-bound P (mol/cm3)
!ELSE        ! deep sea sediments
!   ksPO4 = 3.65D0                                      ! lower sorption rate for deep sea   
                                                       ! Rate constant for kinetic P sorption (1/yr)
!   PO4s = 2.0D-9                                       ! Equilibrium concentration for P sorption (mol/cm3)
!   Minf = 5.2D-15                                      ! asymptotic concentration for Fe-bound P (mol/cm3)
!END IF

!---other parameters
!logdir = "CGEM_log/"

RETURN

END SUBROUTINE default_OMENSED_params

!========================================================================

END MODULE default_OMENSED
