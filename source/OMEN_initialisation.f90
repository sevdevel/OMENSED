MODULE OMEN_initialisation
!************************************************************************
!
! *OMEN_initialisation* routines to initialise the OMENSED model
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Last update - 10 Nov 2021  @(OMENSED)OMEN_initialisation.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - CGEM_model
!
! Internal calls -
!
! Module calls -
!
!************************************************************************

USE OMENpars
USE OMENvars

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE initialise_OMEN_params(WaterTemp, O2_BW, sed_w)
!************************************************************************
!
! *initialise_OMEN_params* Initialise variable parameters for the OMENSED model
!
! Author - Sebastiaan van de Velde
!
! Version - @(OMENSED)initialise_OMEN_params.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - CGEM_model
!
!************************************************************************
!
        
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: WaterTemp, O2_BW, sed_w
!
!*Local variables
!
REAL(8) :: loc_k_apparent
REAL(8) :: loc_total_POC_flux   ! total POC flux at SWI (POC1 + POC2) [mol/(cm^2 yr)]
INTEGER :: i
REAL(8), DIMENSION(OS_vertical_grid) :: OS_z_auxvector

! initialise temperature dependent parameters

DO21=(qdispO2+adispO2*WaterTemp)*dispFactor+Dbio      ! O2 diffusion coefficient in bioturbated layer (cm2/yr)
DO22=(qdispO2+adispO2*WaterTemp)*dispFactor           ! O2 diffusion coefficient in non-bioturbated layer (cm2/yr)
DN1=(qdispNO3+adispNO3*WaterTemp)*dispFactor+Dbio
DN2=(qdispNO3+adispNO3*WaterTemp)*dispFactor
DNH41=((qdispNH4+adispNH4*WaterTemp)*dispFactor+Dbio)/(1.0D0+KNH4)
DNH42=((qdispNH4+adispNH4*WaterTemp)*dispFactor)/(1.0D0+KNH4)
DSO41=(qdispSO4+adispSO4*WaterTemp)*dispFactor+Dbio    ! SO4 diffusion coefficient in bioturbated layer (cm2/yr)
DSO42=(qdispSO4+adispSO4*WaterTemp)*dispFactor         ! SO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
DH2S1=(qdispH2S+adispH2S*WaterTemp)*dispFactor+Dbio
DH2S2=(qdispH2S+adispH2S*WaterTemp)*dispFactor
DPO41=((qdispPO4+adispPO4*WaterTemp)*dispFactor+Dbio)  ! PO4 diffusion coefficient in bioturbated layer (cm2/yr)
DPO42=((qdispPO4+adispPO4*WaterTemp)*dispFactor);      ! PO4 diffusion coefficient in non-bioturbated layer (cm2/yr)
DDIC1=(qdispDIC+adispDIC*WaterTemp)*dispFactor+Dbio    ! DIC diffusion coefficient in bioturbated layer (cm2/yr)
DDIC2=(qdispDIC+adispDIC*WaterTemp)*dispFactor         ! DIC diffusion coefficient in non-bioturbated layer (cm2/yr)
DALK1=(qdispALK+adispALK*WaterTemp)*dispFactor+Dbio;   ! ALK diffusion coefficient in bioturbated layer (cm2/yr)
DALK2=(qdispALK+adispALK*WaterTemp)*dispFactor;        ! ALK diffusion coefficient in non-bioturbated layer (cm2/yr)

IF (O2_BW.LE.OS_BW_O2_anoxia) THEN
    ! decrease bioturbation depth
    zbio = 0.01D0
END IF
         
! initialise depth_vector

DO i=0, (OS_vertical_grid-1)

   OS_z_auxvector(i+1) = 0.1D0*(1.037063D0)**i

END DO
      
OS_z_vector(1) = (OS_z_auxvector(1)/2.D0)

DO i=2, OS_vertical_grid

   OS_z_vector(i) = SUM(OS_z_auxvector(1:(i-1))) + (OS_z_auxvector(i)/2.D0)
   !print*, 'z', i, '=', OS_z_vector(i)
END DO
      
RETURN

END SUBROUTINE initialise_OMEN_params

!========================================================================

END MODULE OMEN_initialisation
