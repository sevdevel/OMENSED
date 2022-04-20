
program MAIN
!************************************************************************
!
! *CGEM* Main model program for the CGEM model
!
! Author - Sebastiaan van de Velde (and original CGEM authors)
!
! Version - @(CGEM)CGEM_model.F90  V0.1
!
! Description -
!
! Calling program -
!
!************************************************************************
!
USE OMENSED_module ! OMENSED module for CGEM
USE OMENvars
USE OMENpars
USE OMEN_readwriteroutines
USE OMENids
USE default_OMENSED

IMPLICIT NONE

!*Local variables

INTEGER :: start_t, end_t, i, output_counter
REAL(8) :: t_real, dtyr
CHARACTER(*), PARAMETER :: logname = "DEV_log.log"
!REAL(8), DIMENSION(8) :: scalefactor=(/ 1.D-3, 5.D-2, 1.D-1, 1.D0, 5.D0, 1.D1, 5.D1, 1.D2/)
REAL(8) :: sedvel = 2.0D-1

output_counter = 0

! Open logfile
OPEN(UNIT=1000, FILE = TRIM(TRIM(logdir)//logname), STATUS = 'REPLACE', ACTION = 'READWRITE')

! Get the start time of the simulation
!start_t = read_clock()

! Allocate OMENSED specific arrays (needs to be moved to CGEM subroutine)
CALL default_OMENSED_params
CALL Allocate_OMEN_Arrays
     
! Write output to file
print_results=.FALSE.
    
! Start the time-loop
!WRITE(1000,FMT=*) "Start time-loop"
DO nt=0, MAXT, DELTI

    ! Print the timestep (and equivalent seconds)
    t_real = DBLE(nt)
    !t_real = t_real/(24.*3600.)
    
    print*, "t_real = ", t_real, "days"
    !WRITE(1000,FMT=*) "t_real = ", t_real, "days"
IF (nt.LE.1) THEN
   
   zbio = 1.0D0                                          ! bioturbation depth (cm)
   Dbio = 0.02D0  
   OS_BW_conds(1,iarr_OS_O2)  = 210.D-9          ! mol cm-3
   OS_BW_conds(1,iarr_OS_NO3) = 9.6D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_NH4) = 0.4D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_SO4) = 28.D-6           ! mol cm-3
   OS_BW_conds(1,iarr_OS_H2S) = 0.D0             ! mol cm-3
   OS_BW_conds(1,iarr_OS_PO4) = 0.0D-9            ! mol cm-3
   OS_BW_conds(1,iarr_OS_DIC) = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_AT)  = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_POC) = 10.D0*1.D-3*1.D-4*365.25D0
   
   OS_BW_conds(1,iarr_OS_bsi) = 10.D0*1.D-3*1.D-4*365.25D0

ELSEIF (nt.LE.2) THEN
   
   zbio = 1.0D-2                                          ! bioturbation depth (cm)
   Dbio = 0.02D0  
   OS_BW_conds(1,iarr_OS_O2)  = 10.D-9          ! mol cm-3
   OS_BW_conds(1,iarr_OS_NO3) = 25.0D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_NH4) = 0.0D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_SO4) = 28.D-6           ! mol cm-3
   OS_BW_conds(1,iarr_OS_H2S) = 0.D0             ! mol cm-3
   OS_BW_conds(1,iarr_OS_PO4) = 50.0D-9            ! mol cm-3
   OS_BW_conds(1,iarr_OS_DIC) = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_AT)  = 2480.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_POC) = 10.D0*1.D-3*1.D-4*365.25D0
   OS_BW_conds(1,iarr_OS_bsi) = 10.D0*1.D-3*1.D-4*365.25D0

ELSEIF (nt.LE.3) THEN

   zbio = 1.0D+1                                          ! bioturbation depth (cm)
   Dbio = 0.17D0  
   OS_BW_conds(1,iarr_OS_O2)  = 250.D-9          ! mol cm-3
   OS_BW_conds(1,iarr_OS_NO3) = 25.0D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_NH4) = 0.6D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_SO4) = 28.D-6           ! mol cm-3
   OS_BW_conds(1,iarr_OS_H2S) = 0.D0             ! mol cm-3
   OS_BW_conds(1,iarr_OS_PO4) = 0.0D-9            ! mol cm-3
   OS_BW_conds(1,iarr_OS_DIC) = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_AT)  = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_POC) = 10.D0*1.D-3*1.D-4*365.25D0
   OS_BW_conds(1,iarr_OS_bsi) = 10.D0*1.D-3*1.D-4*365.25D0

ELSEIF (nt.LE.4) THEN

   zbio = 4.2D0                                          ! bioturbation depth (cm)
   Dbio = 0.18D0  
   OS_BW_conds(1,iarr_OS_O2)  = 243.D-9          ! mol cm-3
   OS_BW_conds(1,iarr_OS_NO3) = 30.1D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_NH4) = 0.22D-9           ! mol cm-3
   OS_BW_conds(1,iarr_OS_SO4) = 28.D-6           ! mol cm-3
   OS_BW_conds(1,iarr_OS_H2S) = 0.D0             ! mol cm-3
   OS_BW_conds(1,iarr_OS_PO4) = 0.0D-9            ! mol cm-3
   OS_BW_conds(1,iarr_OS_DIC) = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_AT)  = 2400.D-9         ! mol cm-3
   OS_BW_conds(1,iarr_OS_POC) = 10.D0*1.D-3*1.D-4*365.25D0
   OS_BW_conds(1,iarr_OS_bsi) = 10.D0*1.D-3*1.D-4*365.25D0
   
END IF
!Call OMENSED

CALL OMENSED(SedVel=sedvel, dum_sfcsumocn=OS_BW_conds(1,:), &
           & dum_new_swifluxes=OS_SWI_fluxes(1,:), dum_OMENSED_BC=OS_boundaryconds(1,:), &
           & dum_diag_profile=OS_diag_profile, dum_z_vector=OS_z_vector, &
           & dum_RCM_approx=OS_RCM_array, dum_bsi_approx=OS_bsi_array, &
           & dum_POC_conc_swi=OS_POC_conc_swi, dum_bsi_conc_swi=OS_bsi_conc_swi) 

  
CALL write_OMENSED_output(filename="output.dat",iounit=2000,fluxarray=OS_SWI_fluxes, &
                        & BCarray=OS_boundaryconds,delxi=1,output_counter=nt)

CALL write_OMENSED_diagprof_output(filename="DiagProfoutput.dat",iounit=2001,diag_profile=OS_diag_profile(:,:),&
                                 & z_vector=OS_z_vector,delxi=1,output_counter=nt)

END DO

! Close logfile

CALL Deallocate_OMEN_Arrays
CLOSE(1000)

end program
