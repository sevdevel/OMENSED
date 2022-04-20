MODULE OMEN_readwriteroutines
!************************************************************************
!
! *OMEN_readwriteroutines* routines to read and write to files for the OMENSED model
!
! Author - Sebastiaan van de Velde
!
! Last update - 10 Nov 2021  @(OMENSED)OMEN_readwriteroutines.f90  V0.1
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
USE OMENids

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE write_OMENSED_output(filename,iounit,fluxarray,BCarray,delxi,output_counter)
!************************************************************************
!
! *write_output* write hydrodynamic output to a file
!
! Author - Sebastiaan van de Velde
!
! Version - @(OMENSED)OMEN_readwriteroutines.f90  V0.1
!
! Description -
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: filename
INTEGER, INTENT(IN) :: iounit, output_counter, delxi
REAL(8), DIMENSION(1,MaxOMENSWIArids), INTENT(IN) :: fluxarray
REAL(8), DIMENSION(1,MaxOMENbcArids), INTENT(IN) :: BCarray

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iounit*    INTEGER File unit
!*output_counter* INTEGER counts how many time output has been written
!*filename*  CHAR    File name
!*array*     REAL    array value to be written to the output file
!*new*       LOGICAL TRUE if first time this file is written
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i
LOGICAL :: Iquestionmyexistence
CHARACTER(40) :: filenameplace

filenameplace = filename!TRIM(outdir)//filename

! Check whether file already exists
INQUIRE(FILE=TRIM(filenameplace), EXIST=Iquestionmyexistence)

IF (Iquestionmyexistence) THEN

    ! If file exists, check if it is because we are already running the model
    IF (output_counter.EQ.0) THEN
        ! If this is first time output has been written during this model run, remove the old file
        OPEN(UNIT=iounit, FILE = TRIM(filenameplace), STATUS = 'REPLACE', ACTION = 'READWRITE')
    ELSE
        ! If this is not the first time output has been written during this model run, add to the existing file
        OPEN(UNIT=iounit, FILE = TRIM(filenameplace), STATUS = 'OLD', ACTION = 'READWRITE', POSITION="append")
    END IF
ELSE
    ! If file does not exist, make one
    OPEN(UNIT=iounit, FILE = TRIM(filenameplace), STATUS = 'NEW', ACTION = 'READWRITE')
END IF

! If this is the first time output is written, add the fluxnames
IF (output_counter.EQ.0) THEN
    WRITE(iounit,FMT=*) 'nt', ' ', 'w', ' ', 'POCpres', ' ', 'zox', ' ', 'zno3', ' ', 'zso4', ' ', &
                      & 'flux_POC', ' ', 'flux_O2', ' ', 'flux_NO3', ' ', &
                      & 'flux_NH4', ' ', 'flux_SO4', ' ', 'flux_H2S', ' ', 'flux_PO4', ' ', 'flux_DIC', ' ', 'flux_AT'
END IF

! Write data from array to file
WRITE(iounit,FMT=*) nt, BCarray(1,iarr_OS_w), BCarray(1,iarr_OS_POC_pres_frac), BCarray(1,iarr_OS_zox), &
                      & BCarray(1,iarr_OS_zno3), BCarray(1,iarr_OS_zso4),  &
                      & fluxarray(1,iarr_OS_POC),  &
                      & fluxarray(1,iarr_OS_O2), fluxarray(1,iarr_OS_NO3), fluxarray(1,iarr_OS_NH4), &
                      & fluxarray(1,iarr_OS_SO4), fluxarray(1,iarr_OS_H2S), fluxarray(1,iarr_OS_PO4), &
                      & fluxarray(1,iarr_OS_DIC), fluxarray(1,iarr_OS_AT)
! Close file
CLOSE(iounit)

END SUBROUTINE write_OMENSED_output

!========================================================================

SUBROUTINE write_OMENSED_diagprof_output(filename,iounit,diag_profile,z_vector,delxi,output_counter)
!************************************************************************
!
! *write_output* write hydrodynamic output to a file
!
! Author - Sebastiaan van de Velde
!
! Version - @(OMENSED)OMEN_readwriteroutines.f90  V0.1
!
! Description -
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: filename
INTEGER, INTENT(IN) :: iounit, output_counter, delxi
REAL(8), DIMENSION(OS_vertical_grid), INTENT(IN) :: z_vector
REAL(8), DIMENSION(OS_vertical_grid,MaxOMENSWIArids), INTENT(IN) :: diag_profile

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iounit*    INTEGER File unit
!*output_counter* INTEGER counts how many time output has been written
!*filename*  CHAR    File name
!*array*     REAL    array value to be written to the output file
!*new*       LOGICAL TRUE if first time this file is written
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j
LOGICAL :: Iquestionmyexistence
CHARACTER(40) :: filenameplace

filenameplace = filename!TRIM(outdir)//filename

! Check whether file already exists
INQUIRE(FILE=TRIM(filenameplace), EXIST=Iquestionmyexistence)

IF (Iquestionmyexistence) THEN

    ! If file exists, check if it is because we are already running the model
    IF (output_counter.EQ.0) THEN
        ! If this is first time output has been written during this model run, remove the old file
        OPEN(UNIT=iounit, FILE = TRIM(filenameplace), STATUS = 'REPLACE', ACTION = 'READWRITE')
    ELSE
        ! If this is not the first time output has been written during this model run, add to the existing file
        OPEN(UNIT=iounit, FILE = TRIM(filenameplace), STATUS = 'OLD', ACTION = 'READWRITE', POSITION="append")
    END IF
ELSE
    ! If file does not exist, make one
    OPEN(UNIT=iounit, FILE = TRIM(filenameplace), STATUS = 'NEW', ACTION = 'READWRITE')
END IF

! If this is the first time output is written, add the depth values
IF (output_counter.EQ.0) THEN
    WRITE(iounit,FMT=*) 'NA', (z_vector(i), i=1, OS_vertical_grid)    
END IF

! Write data from array to file
DO j=1, MaxOMENSWIArids
 
    WRITE(iounit,FMT=*) nt, (diag_profile(i,j), i=1, OS_vertical_grid)                      

END DO

  
! Close file
CLOSE(iounit)

END SUBROUTINE write_OMENSED_diagprof_output

!========================================================================

END MODULE OMEN_readwriteroutines
