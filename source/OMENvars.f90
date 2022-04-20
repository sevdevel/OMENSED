MODULE OMENvars
!************************************************************************
!
! *CGEMpars* CGEM estuarine model arrays
!
! Author - Sebastiaan van de Velde
!
! Last update - 20 Jan 2021  @(CGEM)CGEMvars.f90  V0.1
!
! Description -
!
!************************************************************************
!

IMPLICIT NONE


!---OMEN variables
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: OS_SWI_fluxes, OS_BW_conds, OS_part_fluxes
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: OS_boundaryconds
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: OS_RCM_array 
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: OS_bsi_array 
!REAL(8), ALLOCATABLE, DIMENSION(:) :: OS_POC_profile
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: OS_diag_profile
REAL(8), ALLOCATABLE, DIMENSION(:) :: OS_z_vector
REAL(8), ALLOCATABLE, DIMENSION(:) :: OS_POC_conc_swi,OS_bsi_conc_swi 

REAL(8), ALLOCATABLE, DIMENSION(:) :: aa1, bb1, aa2, A11, A22
!
! Name        Type    Purpose                                                 Unit
!------------------------------------------------------------------------------
!*C*          REAL    coefficients for tridiagonal matrix                     [-]
!*Z*          REAL    coefficients for tridiagonal matrix                     [-]
!*ZZ*         REAL    Cross-section at reference level                        [m2]
!*H*          REAL    Free cross-section                                      [m2]
!*TH*         REAL    Temporary free cross-section                            [m2]
!*D*          REAL    Total cross-section                                     [m2]
!*Dold*       REAL    Total cross-section                                     [m2]
!*Dold2*      REAL    Total cross-section                                     [m2]
!*U*          REAL    Velocity                                                [m s-1]
!*TU*         REAL    Temporary velocity                                      [m s-1]
!*B*          REAL    Width                                                   [m]
!*Chezy*      REAL    Chézy coefficient                                       [m s-2]
!*Y*          REAL    Convergence                                             [-]
!*E*          REAL    Convergence                                             [-]
!*FRIC*       REAL    Friction coefficient                                    [???]
!*DEPTH*      REAL    Water depth                                             [m]
!*slope*      REAL    Depth at reference level                                [m]
!*Detailed_Chezy* REAL array with the Chezy coefficient
!*Detailed_Mero*  REAL array with the Mero coefficient

SAVE

END MODULE OMENvars
