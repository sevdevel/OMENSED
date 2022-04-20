MODULE OMENids
!************************************************************************
!
! *OMENSEDids* OMENSED model variable ids
!
! Author - Sebastiaan van de Velde
!
! Last update - 7 Jan 2009  @(OMENSED)OMENSEDids.f90  V0.1
!
! Description -
!
! Reference -
!
!************************************************************************
!
!USE OMENpars

IMPLICIT NONE

!---start key id number
INTEGER, PARAMETER, PRIVATE :: n0 = 0, nbc = 0

!
!2. Biogeochemical arrays
!--------------------
!
!--state variables
INTEGER, PARAMETER :: &
& iarr_OS_O2=n0+1,iarr_OS_NO3=n0+2,iarr_OS_NH4=n0+3, iarr_OS_SO4=n0+4,iarr_OS_H2S=n0+5,iarr_OS_PO4=n0+6, &
& iarr_OS_DIC=n0+7,iarr_OS_AT=n0+8, iarr_OS_T=n0+9, iarr_OS_M=n0+10, iarr_OS_POC=n0+11
INTEGER,PARAMETER :: &
& iarr_OS_w=nbc+1, iarr_OS_zox=nbc+2, iarr_OS_zno3=nbc+3, iarr_OS_zSO4=nbc+4, &
& iarr_OS_POC_pres_frac=nbc+5, iarr_OS_POP_pres_frac=nbc+6, iarr_OS_POC_swi=nbc+7, iarr_OS_POC_burial=nbc+8

END MODULE OMENids
