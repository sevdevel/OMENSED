MODULE OMEN_Functions
!************************************************************************
!
! *OMEN_Functions* set of functions used in the OMENSED sediment model
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Last update - 09 Nov 2021  @(OMENSED)OMEN_Functions.f90  V0.1
!
! Description -
!
! Reference -
!
! Calling program - OMENSED_module
!
! Internal calls -
!
! Module calls -
!
!************************************************************************

USE OMENpars
USE OMENvars
USE OMENids

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE RCM_approx(RCMarray,RCM_a,RCM_nu)
!************************************************************************
!
! *RCM_approx* RCM approximation 
!
! Author - Sebastiaan van de Velde
!
! Version - @(OMENSED)RCM_approx.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8),DIMENSION(2,OS_RCM_fracs),INTENT(INOUT) :: RCMarray ! Array of deconvoluted RCM
REAL(8), INTENT(IN) :: RCM_a,RCM_nu
!
!*Local variables
!
REAL(8) :: loc_emin, loc_emax, loc_ne, loc_G0in, loc_G1in, loc_G0, loc_G1
REAL(8),DIMENSION(OS_RCM_fracs) :: loc_k, loc_kk, loc_POCfrac
INTEGER :: i

loc_emin = -15.0D0
loc_emax = -DLOG10(RCM_a)+2.0D0 ! upper k limit for k-bins of multi-G approximation
IF (loc_emin.GE.loc_emax) print*, 'emin >= emax, this cannot be!'

loc_k(1)       = 10.D0**(loc_emin)
loc_kk(1)      = 10.D0**(loc_emin)
loc_POCfrac(1) = gammp(RCM_nu,RCM_a*loc_k(1))                      ! lower gamma function
loc_kk(OS_RCM_fracs)      = 10.D0**loc_emax
loc_k(OS_RCM_fracs)       = 10.D0**loc_emax
loc_POCfrac(OS_RCM_fracs) = GAMMQ(RCM_nu,RCM_a*loc_k(OS_RCM_fracs))! upper gamma function

DO i=2, OS_RCM_fracs-1

    loc_ne        = loc_emin+DBLE((i-1))*(loc_emax-loc_emin)/DBLE((OS_RCM_fracs-1))
    loc_kk(i) = 10**loc_ne 
    
    loc_G0in = RCM_a*loc_kk(i-1)
    loc_G1in = RCM_a*loc_kk(i)

    loc_G0 = GAMMQ(RCM_nu,loc_G0in)
    loc_G1 = GAMMQ(RCM_nu,loc_G1in)
    
    loc_POCfrac(i) = loc_G0 - loc_G1
    loc_k(i)       = loc_kk(i-1) + (loc_kk(i)-loc_kk(i-1))/2.D0
   
END DO 

IF (DABS(SUM(loc_POCfrac)-1.D0).GT.1.D-4) print*, 'SUM(F) .NE . 1'

RCMarray(1,:) = loc_POCfrac
RCMarray(2,:) = loc_k

!Fnonbioi = F.* ( swi.C0*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
!C0i = F.*swi.C0;

RETURN

END SUBROUTINE RCM_approx

!========================================================================

FUNCTION GAMMLN(XX)
!************************************************************************
!
! *GAMMLN* Gamma function 
!
! Author - Donald G. Luttermoser, ETSU/Physics, 2 October 2013
!
! Version - @(OMENSED)RCM_approx.f90  V0.1
!
! Description - USES gcf,gser
!               Returns the incomplete gamma function P(a,x)
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************       
         
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: GAMMLN, XX
!
!*Local variables
!
INTEGER :: J
REAL(8) :: SER, STP, TMP, X, Y, COF(6)
SAVE COF, STP

DATA COF, STP/76.18009172947146D0, -86.50532032941677D0, &
            & 24.01409824083091D0, -1.231739572450155D0, &
            & 0.1208650973866179D-2, -0.5395239384953D-5, &
            & 2.5066282746310005D0/

X = XX
Y = X
TMP = X + 5.5D0
TMP = (X + 0.5D0) * DLOG(TMP) - TMP
SER = 1.000000000190015D0

DO 11 J = 1, 6
    Y = Y + 1.D0
    SER = SER + COF(J)/Y
11   CONTINUE

GAMMLN = TMP + DLOG(STP * SER / X)

RETURN

END FUNCTION

!========================================================================
      
FUNCTION GAMMP(a,x)
!************************************************************************
!
! *GAMMP* Gamma function 
!
! Author - Donald G. Luttermoser, ETSU/Physics, 2 October 2013
!
! Version - @(OMENSED)RCM_approx.f90  V0.1
!
! Description - USES gcf,gser
!               Returns the incomplete gamma function P(a,x)
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!      
REAL(8) :: a,x,GAMMP
!
!*Local variables
! 
REAL(8) :: gammcf,gamser,gln

!IF(x.LT.0.D0 .OR. a.LE.0.D0) !pause bad arguments in GAMMP 
!print*,'put',a,x
IF(x.LT.(a+1.D0)) THEN !Use the series representation.
   CALL gser(gamser,a,x,gln)
   GAMMP=gamser
ELSE !Use the continued fraction representation
   CALL gcf(gammcf,a,x,gln)
   GAMMP=1.D0-gammcf ! and take its complement
END IF

RETURN 

END FUNCTION

!========================================================================

FUNCTION GAMMQ(A, X)
!************************************************************************
!
! *GAMMQ* Gamma function 
!
! Author - Donald G. Luttermoser, ETSU/Physics, 2 October 2013
!
! Version - @(OMENSED)RCM_approx.f90  V0.1
!
! Description - USES gcf,gser
!               Returns the incomplete gamma function P(a,x)
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: A, GAMMQ, X
!
!*Local variables
!
REAL(8) :: GAMMCF, GAMSER, GLN

! Check passed parameters.

IF ((X.LT.0.D0).OR.(A.LE.0.D0)) THEN
   print*, 'Bad arguments in GAMMQ.'
   GAMMQ = 0.D0
   RETURN
END IF

IF (X.LT.(A+1.D0)) THEN
   CALL GSER(GAMSER, A, X, GLN)
   GAMMQ = 1.D0 - GAMSER
ELSE
   CALL GCF(GAMMCF, A, X, GLN)
   GAMMQ = GAMMCF
ENDIF

RETURN

END FUNCTION

!========================================================================

SUBROUTINE GCF(GAMMCF, A, X, GLN)
!************************************************************************
!
! *GCF* Gamma function 
!
! Author - Donald G. Luttermoser, ETSU/Physics, 2 October 2013
!
! Version - @(OMENSED)RCM_approx.f90  V0.1
!
! Description - USES gcf,gser
!               Returns the incomplete gamma function P(a,x)
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: A, GAMMCF, GLN, X
!
!*Local variables
!
INTEGER :: ITMAX
REAL(8) :: EPS, FPMIN
PARAMETER (ITMAX=100, EPS=3.D-7, FPMIN=1.D-30)
INTEGER :: I
REAL(8) :: AN, B, C, D, DEL, H!, GAMMLN

GLN = GAMMLN(A)
B = X + 1.D0 - A
C = 1.D0 / FPMIN
D = 1.D0 / B
H = D
DO 11 I = 1, ITMAX
   AN = -I * (I - A)
   B = B + 2.D0
   D = AN*D + B
   IF (DABS(D).LT.FPMIN) D = FPMIN
   C = B + AN/C
   IF (DABS(C).LT.FPMIN) C = FPMIN
   D = 1.D0 / D
   DEL = D * C
   H = H * DEL
   IF (DABS(DEL-1.D0).LT.EPS) GOTO 1
 11   CONTINUE

print*, 'A too large, ITMAX too small in GCF.'
1   GAMMCF = DEXP(-X + A*DLOG(X) - GLN) * H

RETURN

END SUBROUTINE GCF

!========================================================================

SUBROUTINE GSER(GAMSER, A, X, GLN)
!************************************************************************
!
! *GSER* Gamma function 
!
! Author - Donald G. Luttermoser, ETSU/Physics, 2 October 2013
!
! Version - @(OMENSED)RCM_approx.f90  V0.1
!
! Description - USES gcf,gser
!               Returns the incomplete gamma function P(a,x)
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
! 
REAL(8) :: A, GAMSER, GLN, X
!
!*Local variables
!
REAL(8) :: EPS
INTEGER :: ITMAX
PARAMETER (ITMAX=100, EPS=3.D-7)
INTEGER :: N
REAL(8) :: AP, DEL, SUM!, GAMMLN

! Call the GAMMLN function.

!print*,'From GSER',A,X
GLN = GAMMLN(A)
IF (X.LE.0.D0) THEN
   IF (X.LT.0.D0) print*, 'X < 0 in GSER.'
   GAMSER = 0.D0
   RETURN
END IF

AP = A
SUM = 1.D0/A
DEL = SUM

DO 11 N = 1, ITMAX
   AP = AP + 1.D0
   DEL = DEL * X / AP
   SUM = SUM + DEL
   IF (DABS(DEL).LT.DABS(SUM)*EPS) GOTO 1
11   CONTINUE

print*, 'A too large, ITMAX too small in GSER.'

1   GAMSER = SUM * DEXP(-X + A*DLOG(X) - GLN)

RETURN
      
END SUBROUTINE GSER
      
!========================================================================
      
SUBROUTINE Conservation_Check(fluxarray,bcarray)
!************************************************************************
!
! *Conservation_Check* check conservation for elements and electrons 
!
! Author - Sebastiaan van de Velde
!
! Version - @(OMENSED)Conservation_Check.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - OMENSED_module
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8),DIMENSION(MaxOMENSWIArids),INTENT(IN) :: fluxarray ! SWI fluxes
REAL(8),DIMENSION(MaxOMENbcArids),INTENT(IN) :: bcarray   ! other vars
!
!*Local variables
!
REAL(8) :: tol  ! imbalance tolerance
REAl(8) :: loc_Ccheck, loc_Scheck, loc_Pcheck, loc_Ncheck, loc_echeck
REAl(8) :: loc_POC_remin, loc_POP_remin

tol           = 1.D-20
loc_POC_remin = -(1.D0-bcarray(iarr_OS_POC_pres_frac))*fluxarray(iarr_OS_POC)
loc_POP_remin = -(1.D0-bcarray(iarr_OS_POP_pres_frac))*Z_P/X_C*fluxarray(iarr_OS_POC)

! Carbon check
loc_Ccheck = loc_POC_remin + &                               ! remineralised POC flux
           & fluxarray(iarr_OS_DIC)                          ! DIC efflux 

! Sulfur check          
loc_Scheck = fluxarray(iarr_OS_SO4) + fluxarray(iarr_OS_H2S) ! S balance 

! Phosphorus check
loc_Pcheck = loc_POP_remin + &                           ! remineralised POP flux
           & fluxarray(iarr_OS_PO4)                          ! PO4 efflux 

! Nitrogen check
loc_Ncheck = Y_N/X_C*loc_POC_remin + &  ! remineralised PON flux
           & fluxarray(iarr_OS_NH4) + fluxarray(iarr_OS_NO3)      ! NH4 efflux 
         
! electron check
loc_echeck = 4.0D0*loc_POC_remin                                  + & ! electrons from C-CH2O  to C-CO2 (-4)
           & 4.0D0*(Y_N/X_C*loc_POC_remin+fluxarray(iarr_OS_NH4)) - & ! electrons from N-NH4   to N-NO3 (-4)
           & 4.0D0*fluxarray(iarr_OS_O2)                          - & ! electrons from O-O2    to O-H2O (+4) 
           & 5.0D0*fluxarray(iarr_OS_NO3)                         - & ! electrons from N-NO3-  to N-N2  (+5) 
           & 8.0D0*fluxarray(iarr_OS_SO4)                             ! electrons from S-SO42- to S-H2S (+8)
           
           
print*, "Ccheck", loc_Ccheck*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
print*, "Scheck", loc_Scheck*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
print*, "Pcheck", loc_Pcheck*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
print*, "Ncheck", loc_Ncheck*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'
print*, "echeck", loc_echeck*1.D3*1.D4/365.25D0, 'mmol m-2 d-1'

RETURN

END SUBROUTINE Conservation_Check

!========================================================================

SUBROUTINE sub_calcfg_l12(z, reac, ktemp, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ls_D1, ls_D2,&
                        & ltype, e, dedz, f, dfdz, g, dgdz, dum_POC_conc_swi, dum_RCM_approx) ! MULTIG
!************************************************************************
!
! *sub_calcfg_l12* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - calculate solution basis functions, for layer which may cross bioturbation boundary
!
!               reac1, reac2        - mol/mol S released per organic carbon C
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!
!               Where e,f,g are generated by matching solutions across bioturbation boundary (if necessary)
!               Solution properties (matching etc) are input in ls
!               On input, ls should contain fields generated by sub_prepfg_l12
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
REAL(8), INTENT(IN) :: z, reac, ktemp
INTEGER, INTENT(IN) :: ltype
REAL(8), INTENT(INOUT) :: e, dedz, f, dfdz, g, dgdz
REAL(8), INTENT(IN) :: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ls_D1, ls_D2
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1

SELECT CASE (ltype)
            
    CASE (1)    ! bioturbated
       CALL sub_calcfg_l1(z=z, reac=reac, Dtemp=ls_D1, ktemp=ktemp, e=e, dedz=dedz, f=f, dfdz=dfdz, &
                        & g=g, dgdz=dgdz, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
    CASE (2)    ! not bioturbated
       CALL sub_calcfg_l2(z=z, reac=reac, Dtemp=ls_D2, ktemp=ktemp, e=e, dedz=dedz, f=f, dfdz=dfdz, &
                        & g=g, dgdz=dgdz)
    CASE (3)    ! crossing boundary
       IF (z.GE.zbio) THEN      ! below bioturbated region
          CALL sub_calcfg_l2(z=z, reac=reac, Dtemp=ls_D2, ktemp=ktemp, e=e, dedz=dedz, f=f, dfdz=dfdz, &
                           & g=g, dgdz=dgdz)
       ELSE    ! above bioturbated region
          CALL sub_calcfg_l1(z=z, reac=reac, Dtemp=ls_D1, ktemp=ktemp, e=e_1, dedz=dedz_1, f=f_1, &
                           & dfdz=dfdz_1, g=g_1, dgdz=dgdz_1, dum_POC_conc_swi=dum_POC_conc_swi, &
                           & dum_RCM_approx=dum_RCM_approx)
          ! Now find 'transformed' basis functions such that in layer 1, O2 = A_2*et + B_2*ft + gt
          ! (ie layer 1 soln written in terms of layer 2 coeffs A_2, B_2)
          CALL sub_xformsoln(e_1, f_1, g_1, dedz_1, dfdz_1, dgdz_1, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, &
                           & e, f, g, dedz, dfdz, dgdz)
       END IF
    CASE DEFAULT
       STOP
END SELECT

END SUBROUTINE sub_calcfg_l12

!========================================================================

SUBROUTINE sub_calcfg_l1(z, reac, Dtemp, ktemp, e, dedz, f, dfdz, g, dgdz, & ! MULTIG
                       & dum_POC_conc_swi, dum_RCM_approx)
!************************************************************************
!
! *sub_calcfg_l1* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Basis functions for solutes, case z <= zbio
!
!               reac1, reac2        - mol./mol S released per organic carbon C
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: z, reac, Dtemp, ktemp                                                          ! in from SUBROUTINE before
REAL(8), INTENT(INOUT) :: e, dedz, f, dfdz, g, dgdz             ! out
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: b1, pfac
REAL(8), DIMENSION(OS_RCM_fracs) :: PhiI, PhiII, PhiIII
REAL(8), DIMENSION(OS_RCM_fracs) :: ea1z, eb1z
REAL(8), DIMENSION(OS_RCM_fracs) :: loc_k
INTEGER :: i

e = 1.0
dedz = 0.0
b1=w/(Dtemp + OS_tolcte)
f=DEXP(z*b1)
dfdz = b1*DEXP(z*b1)

pfac = 1.0D0                    ! in fact, already has (1-por)/por

loc_k(:) = dum_RCM_approx(2,:)

DO i=1, OS_RCM_fracs

   PhiI(i)   = -pfac*loc_k(i)*(reac)*A11(i)/(Dtemp*aa1(i)**2.D0-w*aa1(i)-ktemp + &
             & OS_tolcte)
   PhiII(i)  = pfac*loc_k(i)*(reac)*A11(i)/(Dtemp*bb1(i)**2.D0-w*bb1(i)-ktemp + &
             & OS_tolcte)
   PhiIII(i) = -pfac*loc_k(i)*(reac)*dum_POC_conc_swi(i)/(Dtemp*bb1(i)**2.D0-w*bb1(i)-ktemp &
             & + OS_tolcte)

   ea1z(i) = DEXP(aa1(i)*z)
   eb1z(i) = DEXP(bb1(i)*z)

END DO

g = 0.D0
dgdz = 0.D0

DO i=1, OS_RCM_fracs
   g    = g +(PhiI(i)*ea1z(i) + PhiII(i)*eb1z(i) + PhiIII(i)*eb1z(i))
   dgdz = dgdz + (PhiI(i)*aa1(i)*ea1z(i) + PhiII(i)*bb1(i)*eb1z(i) + PhiIII(i)*bb1(i)*eb1z(i))
END DO

END SUBROUTINE sub_calcfg_l1

!========================================================================

SUBROUTINE sub_calcfg_l2(z, reac, Dtemp, ktemp, e, dedz, f, dfdz, g, dgdz) ! MULTIG
!************************************************************************
!
! *sub_calcfg_l2* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Basis functions for solutes, case z <= zbio
!
!               reac1, reac2        - mol./mol S released per organic carbon C
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: z, reac, Dtemp, ktemp                                                          ! in from SUBROUTINE before
REAL(8), INTENT(INOUT) :: e, dedz, f, dfdz, g, dgdz             ! out
!
!*Local variables
!
REAL(8) :: b2, pfac
REAL(8), DIMENSION(OS_RCM_fracs) :: PhiI, loc_k
INTEGER :: i

e = 1.0D0
dedz = 0.0D0
b2 = w/(Dtemp + OS_tolcte)
f = DEXP(z*b2)
dfdz = b2*DEXP(z*b2)

pfac = 1            !in fact, already has (1-por)/por

loc_k(:) = OS_RCM_array(2,:)

DO i=1, OS_RCM_fracs
   PhiI(i) = -pfac*loc_k(i)*(reac)*A22(i)/(Dtemp*aa2(i)**2.D0-w*aa2(i)-ktemp + &
            & OS_tolcte)
END DO

g = 0.D0
dgdz = 0.D0

DO i=1, OS_RCM_fracs
   g    = g + PhiI(i)*DEXP(aa2(i)*z)
   dgdz = dgdz + PhiI(i)*aa2(i)*DEXP(aa2(i)*z)
END DO

END SUBROUTINE sub_calcfg_l2

!========================================================================

SUBROUTINE sub_prepfg_l12(reac, ktemp, zU, zL, D1, D2, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f, ltype, & ! MULTIG
                        & dum_POC_conc_swi, dum_RCM_approx)
!************************************************************************
!
! *sub_prepfg_l12* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: reac, ktemp, zU, zL, D1, D2
REAL(8), INTENT(INOUT) :: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
INTEGER, INTENT(INOUT) :: ltype
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: e_zbio_l1, dedz_zbio_l1, f_zbio_l1, dfdz_zbio_l1, g_zbio_l1, dgdz_zbio_l1
REAL(8) :: e_zbio_l2, dedz_zbio_l2, f_zbio_l2, dfdz_zbio_l2, g_zbio_l2, dgdz_zbio_l2

IF (zL.LE.zbio) THEN  ! wholly within bioturbated layer
   ltype = 1
ELSEIF (zU .GE. zbio) THEN ! wholly within non-bioturbated layer
   ltype = 2
ELSE             ! crossing boundary - sort out solution matching at zbio
   ltype = 3
   CALL sub_calcfg_l1(z=zbio, reac=reac, Dtemp=D1, ktemp=ktemp, e=e_zbio_l1, dedz=dedz_zbio_l1, &
                    & f=f_zbio_l1, dfdz=dfdz_zbio_l1, g=g_zbio_l1, dgdz=dgdz_zbio_l1, &
                    & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                    
   CALL sub_calcfg_l2(z=zbio, reac=reac, Dtemp=D2, ktemp=ktemp, e=e_zbio_l2, dedz=dedz_zbio_l2, &
                    & f=f_zbio_l2, dfdz=dfdz_zbio_l2, g=g_zbio_l2, dgdz=dgdz_zbio_l2)

   ! match solutions at zbio - continuous concentration and flux
   CALL sub_matchsoln(e_zbio_l1, f_zbio_l1, g_zbio_l1, D1*dedz_zbio_l1, D1*dfdz_zbio_l1, &
                    & D1*dgdz_zbio_l1, &
                    & e_zbio_l2, f_zbio_l2, g_zbio_l2, D2*dedz_zbio_l2, D2*dfdz_zbio_l2, D2*dgdz_zbio_l2, &
                    & 0.0D0, 0.0D0, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)
END IF

END SUBROUTINE sub_prepfg_l12

!========================================================================

SUBROUTINE sub_calcfg_l12_PO4_M(z, reac1P, dum_ktempP, dum_QtempP, &
       & dum_D1P, dum_D2P, dum_alphaP, dum_mat_C, dum_vec_D, dum_ltype, &
       & dum_ktempM, dum_QtempM, dum_D1M, dum_D2M, dum_alphaM, e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, &
       & p_P, dpdz_P, q_P, dqdz_P, e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M, &
       & dum_POC_conc_swi, dum_RCM_approx)
!************************************************************************
!
! *sub_calcfg_l12_PO4_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Utility subroutines for PO4/Fe-bound P
!               calculate solution basis functions, for layer which may cross bioturbation boundary
!
!               reac1, reac2        - mol/mol S released per organic carbon C
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!
!               Where e,f,g are generated by matching solutions across bioturbation boundary (if necessary)
!               Solution properties (matching etc) are input in ls
!               On input, ls should contain fields generated by prepfg_l12
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: z, reac1P, dum_ktempP, dum_QtempP, dum_alphaP, dum_D1P, dum_D2P
REAL(8), INTENT(IN) :: dum_ktempM, dum_QtempM, dum_alphaM, dum_D1M, dum_D2M
INTEGER, INTENT(IN) :: dum_ltype
REAL(8), DIMENSION(4,4), INTENT(IN) :: dum_mat_C
REAL(8), DIMENSION(1:4), INTENT(IN) ::  dum_vec_D
REAL(8), INTENT(INOUT) :: e_P, dedz_P, f_P, dfdz_P, g_P, dgdz_P, p_P, dpdz_P, q_P, dqdz_P ! ODE solutions (E, F, P, Q) and the particulat integral (G) and their derivatives
REAL(8), INTENT(INOUT) :: e_M, dedz_M, f_M, dfdz_M, g_M, dgdz_M, p_M, dpdz_M, q_M, dqdz_M
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: loc_a1_P, loc_b1_P, loc_a1_M, loc_b1_M, loc_a2_M, loc_a2_P, loc_b2_P
REAL(8), DIMENSION(3,OS_RCM_fracs) :: loc_Phi1_P
REAL(8), DIMENSION(OS_RCM_fracs) :: loc_Phi2_P
REAL(8) :: loc_e_P_1, loc_dedz_P_1, loc_f_P_1, loc_dfdz_P_1, loc_g_P_1, loc_dgdz_P_1
REAL(8) :: loc_p_P_1, loc_dpdz_P_1, loc_q_P_1, loc_dqdz_P_1
REAL(8) :: loc_e_M_1, loc_dedz_M_1, loc_f_M_1, loc_dfdz_M_1, loc_g_M_1, loc_dgdz_M_1
REAL(8) :: loc_p_M_1, loc_dpdz_M_1, loc_q_M_1, loc_dqdz_M_1
REAL(8), DIMENSION(1:4) :: loc_EFPQ_P, loc_dEFPQdz_P, loc_EFPQ_M, loc_dEFPQdz_M ! save the ODE solutions in vectors to make calculation easier (DH?: however, is this faster?)
REAL(8), DIMENSION(1:4) :: loc_EFPQ_P_t, loc_dEFPQdz_P_t, loc_EFPQ_M_t, loc_dEFPQdz_M_t ! the transformed ODE solutions coming from sub_xformsoln_PO4_M

loc_Phi1_P = 0.D0!(/ 0, 0, 0, 0, 0, 0 /)
loc_Phi2_P = 0.D0!(/ 0, 0, 0, 0, 0, 0 /)

SELECT CASE(dum_ltype)
     CASE(1)    ! bioturbated
       IF (dum_alphaP.EQ.0) THEN   ! oxic layer -> call PO4 first
           CALL sub_calcfg_l1_PO4(z=z, reac=reac1P, Dtemp=dum_D1P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                                & a1_M=0.0D0, b1_M=0.0D0, alpha=dum_alphaP, e=e_P, dedz=dedz_P, f=f_P, &
                                & dfdz=dfdz_P, g=g_P, dgdz=dgdz_P, p=p_P, dpdz=dpdz_P, q=q_P, dqdz=dqdz_P, &
                                & a1=loc_a1_P, b1=loc_b1_P, Phi1=loc_Phi1_P, &
                                & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
           CALL sub_calcfg_l1_M(z=z, Dtemp=dum_D1M, ktemp=dum_ktempM, Qtemp=dum_QtempM, a1_P=loc_a1_P, b1_P=loc_b1_P, &
                              & Phi1_P=loc_Phi1_P, alpha=dum_alphaM, e=e_M, dedz=dedz_M, &
                              & f=f_M, dfdz=dfdz_M, g=g_M, dgdz=dgdz_M, p=p_M, &
                              & dpdz=dpdz_M, q=q_M, dqdz=dqdz_M, c1=loc_a1_M, d1=loc_b1_M)
       ELSE        ! anoxic layer -> call M first
           CALL sub_calcfg_l1_M(z=z, Dtemp=dum_D1M, ktemp=dum_ktempM, Qtemp=dum_QtempM, a1_P=0.0D0, b1_P=0.0D0, &
                                & Phi1_P=loc_Phi1_P, alpha=dum_alphaM, e=e_M, dedz=dedz_M, &
                                & f=f_M, dfdz=dfdz_M, g=g_M, dgdz=dgdz_M, p=p_M, &
                                & dpdz=dpdz_M, q=q_M, dqdz=dqdz_M, c1=loc_a1_M, d1=loc_b1_M)
           CALL sub_calcfg_l1_PO4(z=z, reac=reac1P, Dtemp=dum_D1P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                                & a1_M=loc_a1_M, b1_M=loc_b1_M, alpha=dum_alphaP, e=e_P, dedz=dedz_P, f=f_P, &
                                & dfdz=dfdz_P, g=g_P, dgdz=dgdz_P, p=p_P, dpdz=dpdz_P, q=q_P, dqdz=dqdz_P, &
                                & a1=loc_a1_P, b1=loc_b1_P, Phi1=loc_Phi1_P, &
                                & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
       END IF

     CASE(2)    ! not bioturbated
       IF (dum_alphaP.EQ.0) THEN   ! oxic layer -> call PO4 first
          CALL sub_calcfg_l2_PO4(z=z, reac=reac1P, Dtemp=dum_D2P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                               & a2_M=0.0D0, alpha=dum_alphaP, e=e_P, dedz=dedz_P, f=f_P, dfdz=dfdz_P, &
                               & g=g_P, dgdz=dgdz_P, p=p_P, dpdz=dpdz_P, q=q_P, dqdz=dqdz_P, a2=loc_a2_P, &
                               & b2=loc_b2_P, Phi2=loc_Phi2_P, dum_RCM_approx=dum_RCM_approx)
          CALL sub_calcfg_l2_M(z=z, ktemp=dum_ktempM, Qtemp=dum_QtempM, a2_P=loc_a2_P, b2_P=loc_b2_P, Phi2_P=loc_Phi2_P, &
                                & alpha=dum_alphaM, e=e_M, dedz=dedz_M, f=f_M, dfdz=dfdz_M, g=g_M, dgdz=dgdz_M, p=p_M, &
                                & dpdz=dpdz_M, q=q_M, dqdz=dqdz_M, c2=loc_a2_M)
       ELSE    ! anoxic layer -> call M first
          CALL sub_calcfg_l2_M(z=z, ktemp=dum_ktempM, Qtemp=dum_QtempM, a2_P=0.0D0, b2_P=0.0D0, Phi2_P=loc_Phi2_P, &
                                & alpha=dum_alphaM, e=e_M, dedz=dedz_M, f=f_M, dfdz=dfdz_M, g=g_M, dgdz=dgdz_M, p=p_M, &
                                & dpdz=dpdz_M, q=q_M, dqdz=dqdz_M, c2=loc_a2_M)
          CALL sub_calcfg_l2_PO4(z=z, reac=reac1P, Dtemp=dum_D2P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                               & a2_M=loc_a2_M, alpha=dum_alphaP, e=e_P, dedz=dedz_P, f=f_P, dfdz=dfdz_P, &
                               & g=g_P, dgdz=dgdz_P, p=p_P, dpdz=dpdz_P, q=q_P, dqdz=dqdz_P, a2=loc_a2_P, &
                               & b2=loc_b2_P, Phi2=loc_Phi2_P, dum_RCM_approx=dum_RCM_approx)
       END IF

     CASE(3)    ! crossing boundary
       IF (z.GT.zbio) THEN      ! not bioturbated region
          IF (dum_alphaP.EQ.0) THEN   ! oxic layer -> call PO4 first NOTE: BUT DECIDE VIA ALPHA_M NOT WITH <= ZOX!!! DOESN't WORK FOR BOUNDARY ZOX
             CALL sub_calcfg_l2_PO4(z=z, reac=reac1P, Dtemp=dum_D2P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                                  & a2_M=0.0D0, alpha=dum_alphaP, e=e_P, dedz=dedz_P, f=f_P, dfdz=dfdz_P, &
                                  & g=g_P, dgdz=dgdz_P, p=p_P, dpdz=dpdz_P, q=q_P, dqdz=dqdz_P, a2=loc_a2_P, &
                                  & b2=loc_b2_P, Phi2=loc_Phi2_P, dum_RCM_approx=dum_RCM_approx)
             CALL sub_calcfg_l2_M(z=z, ktemp=dum_ktempM, Qtemp=dum_QtempM, a2_P=loc_a2_P, b2_P=loc_b2_P, Phi2_P=loc_Phi2_P, &
                                & alpha=dum_alphaM, e=e_M, dedz=dedz_M, f=f_M, dfdz=dfdz_M, g=g_M, dgdz=dgdz_M, p=p_M, &
                                & dpdz=dpdz_M, q=q_M, dqdz=dqdz_M, c2=loc_a2_M)
          ELSE ! anoxic layer -> call M first
             CALL sub_calcfg_l2_M(z=z, ktemp=dum_ktempM, Qtemp=dum_QtempM, a2_P=0.0D0, b2_P=0.0D0, Phi2_P=loc_Phi2_P, &
                                & alpha=dum_alphaM, e=e_M, dedz=dedz_M, f=f_M, dfdz=dfdz_M, g=g_M, dgdz=dgdz_M, p=p_M, &
                                & dpdz=dpdz_M, q=q_M, dqdz=dqdz_M, c2=loc_a2_M)
             CALL sub_calcfg_l2_PO4(z=z, reac=reac1P, Dtemp=dum_D2P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                                  & a2_M=loc_a2_M, alpha=dum_alphaP, e=e_P, dedz=dedz_P, f=f_P, dfdz=dfdz_P, &
                                  & g=g_P, dgdz=dgdz_P, p=p_P, dpdz=dpdz_P, q=q_P, dqdz=dqdz_P, a2=loc_a2_P, &
                                  & b2=loc_b2_P, Phi2=loc_Phi2_P, dum_RCM_approx=dum_RCM_approx)

          END IF ! (dum_alphaP==0)
       ELSE    ! bioturbated region z <= zbio
          IF (dum_alphaP.EQ.0) THEN   ! oxic layer -> call PO4 first
             ! CASE 1 & 2: LAYER 1: have 4 int. const.
             CALL sub_calcfg_l1_PO4(z=z, reac=reac1P, Dtemp=dum_D1P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                                  & a1_M=0.0D0, b1_M=0.0D0, alpha=dum_alphaP, e=loc_e_P_1, dedz=loc_dedz_P_1, &
                                  & f=loc_f_P_1, dfdz=loc_dfdz_P_1, g=loc_g_P_1, dgdz=loc_dgdz_P_1, p=loc_p_P_1, &
                                  & dpdz=loc_dpdz_P_1, q=loc_q_P_1, dqdz=loc_dqdz_P_1, a1=loc_a1_P, b1=loc_b1_P, &
                                  & Phi1=loc_Phi1_P,dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
             CALL sub_calcfg_l1_M(z=z, Dtemp=dum_D1M, ktemp=dum_ktempM, Qtemp=dum_QtempM, a1_P=loc_a1_P, b1_P=loc_b1_P, &
                                & Phi1_P=loc_Phi1_P, alpha=dum_alphaM, e=loc_e_M_1, dedz=loc_dedz_M_1, &
                                & f=loc_f_M_1, dfdz=loc_dfdz_M_1, g=loc_g_M_1, dgdz=loc_dgdz_M_1, p=loc_p_M_1, &
                                & dpdz=loc_dpdz_M_1, q=loc_q_M_1, dqdz=loc_dqdz_M_1, c1=loc_a1_M, d1=loc_b1_M)
           ! DH: FOR CASE 2: DON'T HAVE D FROM LAYER 2
           ELSE    ! anoxic layer -> call M first
              ! DH: CASE 1: LAYER 2: have 4 int. const.
              CALL sub_calcfg_l1_M(z=z, Dtemp=dum_D1M, ktemp=dum_ktempM, Qtemp=dum_QtempM, a1_P=0.0D0, b1_P=0.0D0, &
                                 & Phi1_P=loc_Phi1_P, alpha=dum_alphaM, e=loc_e_M_1, dedz=loc_dedz_M_1, &
                                 & f=loc_f_M_1, dfdz=loc_dfdz_M_1, g=loc_g_M_1, dgdz=loc_dgdz_M_1, p=loc_p_M_1, &
                                 & dpdz=loc_dpdz_M_1, q=loc_q_M_1, dqdz=loc_dqdz_M_1, c1=loc_a1_M, d1=loc_b1_M)
              CALL sub_calcfg_l1_PO4(z=z, reac=reac1P, Dtemp=dum_D1P, ktemp=dum_ktempP, Qtemp=dum_QtempP, &
                                 & a1_M=loc_a1_M, b1_M=loc_b1_M, alpha=dum_alphaP, e=loc_e_P_1, dedz=loc_dedz_P_1, &
                                 & f=loc_f_P_1, dfdz=loc_dfdz_P_1, g=loc_g_P_1, dgdz=loc_dgdz_P_1, p=loc_p_P_1, &
                                 & dpdz=loc_dpdz_P_1, q=loc_q_P_1, dqdz=loc_dqdz_P_1, a1=loc_a1_P, b1=loc_b1_P, &
                                 & Phi1=loc_Phi1_P, dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
                                 
           END IF  ! (dum_alphaP==0)

           ! Now find 'transformed' basis functions such that in layer 1,
           ! O2 = A_2*et + B_2*ft + gt  (ie layer 1 solution written in terms of layer 2 coeffs A_2, B_2)

           loc_EFPQ_P = (/ loc_e_P_1, loc_f_P_1, loc_p_P_1, loc_q_P_1 /)
           loc_dEFPQdz_P = (/ loc_dedz_P_1, loc_dfdz_P_1, loc_dpdz_P_1, loc_dqdz_P_1 /)
           loc_EFPQ_M = (/ loc_e_M_1, loc_f_M_1, loc_p_M_1, loc_q_M_1 /)
           loc_dEFPQdz_M = (/loc_dedz_M_1, loc_dfdz_M_1, loc_dpdz_M_1, loc_dqdz_M_1/)

           CALL sub_xformsoln_PO4_M(loc_EFPQ_P, loc_EFPQ_M, loc_dEFPQdz_P, loc_dEFPQdz_M, &
                                  & loc_g_P_1, loc_g_M_1,loc_dgdz_P_1, loc_dgdz_M_1, dum_mat_C, dum_vec_D, &
                                  & loc_EFPQ_P_t, g_P, loc_dEFPQdz_P_t, dgdz_P, loc_EFPQ_M_t, &
                                  & g_M, loc_dEFPQdz_M_t, dgdz_M)

           ! WHEN lTYPE=3 - DEAL WITH ONE VARIABLE SHORT FROM LAYER BELOW
           ! FOR CASE 1 & 2: deal with missing values from layer below
           IF (zox.LE.zbio) THEN   ! CASE 1: no F from layer below - DH: no e_M & f_M as well !? but 0 anyway at the moment
              e_P = loc_EFPQ_P_t(1)
              f_P = loc_EFPQ_P_t(2)
              p_P = loc_EFPQ_P_t(3)
              q_P = loc_q_P_1

              dedz_P = loc_dEFPQdz_P_t(1)
              dfdz_P = loc_dEFPQdz_P_t(2)
              dpdz_P = loc_dEFPQdz_P_t(3)
              dqdz_P = loc_dqdz_P_1

              e_M = loc_EFPQ_M_t(1)
              f_M = loc_EFPQ_M_t(2)
              p_M = loc_EFPQ_M_t(3)
              q_M = loc_q_M_1

              dedz_M = loc_dEFPQdz_M_t(1)
              dfdz_M = loc_dEFPQdz_M_t(2)
              dpdz_M = loc_dEFPQdz_M_t(3)
              dqdz_M = loc_dqdz_M_1

            ELSE            ! CASE 2: no Q from layer below - DH: no p_P as well !?
              e_P = loc_EFPQ_P_t(1)
              f_P = loc_EFPQ_P_t(2)
              p_P = loc_p_P_1            
              q_P = loc_q_P_1

              dedz_P = loc_dEFPQdz_P_t(1)
              dfdz_P = loc_dEFPQdz_P_t(2)
              dpdz_P = loc_dpdz_P_1       
              dqdz_P = loc_dqdz_P_1

              e_M = loc_EFPQ_M_t(1)
              f_M = loc_EFPQ_M_t(2)
              p_M = loc_EFPQ_M_t(3)
              q_M = loc_q_M_1

              dedz_M = loc_dEFPQdz_M_t(1)
              dfdz_M = loc_dEFPQdz_M_t(2)
              dpdz_M = loc_dEFPQdz_M_t(3)
              dqdz_M = loc_dqdz_M_1

            END IF ! (zox .LE. zbio)

       END IF  ! (z > zbio)
       
     CASE DEFAULT
        print*, ' unrecognized ltype in  sub_calcfg_l2_PO4 ', dum_ltype
        STOP
END SELECT

END SUBROUTINE sub_calcfg_l12_PO4_M

!========================================================================

SUBROUTINE sub_calcfg_l1_PO4(z, reac, Dtemp, ktemp, Qtemp, a1_M, b1_M, alpha, e, &
                           & dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a1, b1, Phi1, &
                           & dum_POC_conc_swi, dum_RCM_approx)
!************************************************************************
!
! *sub_calcfg_l1_PO4* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Basis functions for solutes, case z <= zbio
!
!               reac1, reac2        - mol./mol S released per organic carbon C
!               depend1,   depend2 coming from other species
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!               and for dependent species
!               S(z) = A .* e(z) + B .* f(z) + C .* p(z) +  D.* q(z) + g(z)
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: z, reac, Dtemp, ktemp, Qtemp, a1_M, b1_M, alpha              ! in from SUBROUTINE before
REAL(8), INTENT(INOUT) ::  e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a1, b1 ! out
!REAL(8), DIMENSION(1:6), INTENT(INOUT) :: Phi1                               ! out
REAL(8), DIMENSION(3,OS_RCM_fracs), INTENT(INOUT) :: Phi1                               ! out
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: pfac
REAL(8) :: ea1z, eb1z
INTEGER :: i, j
REAL(8), DIMENSION(OS_RCM_fracs) :: loc_k

a1 = (w-DSQRT(w**2.0D0+4.0D0*Dtemp*ktemp))/(2.0D0*Dtemp + OS_tolcte)
e = DEXP(z*a1)
dedz = a1*DEXP(z*a1)

b1 = (w+DSQRT(w**2.0D0+4.0D0*Dtemp*ktemp))/(2.0D0*Dtemp + OS_tolcte)
f = DEXP(z*b1)
dfdz = b1*DEXP(z*b1)

pfac = 1.D0                    ! in fact, already has (1-por)/por

loc_k(:) = dum_RCM_approx(2,:)

! NOW to OM reaction terms!
! save all Phis in one variable to pass back

DO i=1, OS_RCM_fracs

   Phi1(1,i) = -pfac*loc_k(i)*(reac)*A11(i)/(Dtemp*aa1(i)**2.0D0-w*aa1(i)-ktemp + OS_tolcte)
   
   Phi1(2,i) =  pfac*loc_k(i)*(reac)*A11(i)/(Dtemp*bb1(i)**2.0D0-w*bb1(i)-ktemp + OS_tolcte)
   
   Phi1(3,i) = -pfac*loc_k(i)*(reac)*dum_POC_conc_swi(i)/(Dtemp*bb1(i)**2.0D0-w*bb1(i)-ktemp + OS_tolcte)
   
   ea1z = DEXP(aa1(i)*z)
   eb1z = DEXP(bb1(i)*z)

END DO

g   =0.D0
dgdz=0.D0

DO i=1, OS_RCM_fracs
    g = g + (Phi1(1,i)*ea1z + Phi1(2,i)*eb1z + Phi1(3,i)*eb1z)
END DO

IF (ktemp.NE.0) THEN        !CHECK: actually no need as ktemp always <> 0
   g =  g + Qtemp/(ktemp + OS_tolcte)   ! here problem if ktemp=0
END IF

DO i=1, OS_RCM_fracs
   dgdz = dgdz + (Phi1(1,i)*aa1(i)*ea1z + Phi1(2,i)*bb1(i)*eb1z + Phi1(3,i)*bb1(i)*eb1z)
END DO

IF (alpha.EQ.0.0D0) THEN      ! was z<=res.zox PO4 is independent of M (no info in alpha)
   p = 0.0D0
   dpdz = 0.0D0
   q = 0.0D0
   dqdz = 0.0D0
ELSE                    ! PO4 is dependent on M
   p = -alpha/(Dtemp*a1_M**2.0D0-w*a1_M-ktemp + OS_tolcte)*DEXP(z*a1_M)
   dpdz = a1_M*p
   q = -alpha/(Dtemp*b1_M**2.0D0-w*b1_M-ktemp + OS_tolcte)*DEXP(z*b1_M)
   dqdz = b1_M*q
END IF

END SUBROUTINE sub_calcfg_l1_PO4

!========================================================================

SUBROUTINE sub_calcfg_l2_PO4(z, reac, Dtemp, ktemp, Qtemp, a2_M, alpha, &
                           & e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a2, b2, Phi2, &
                           & dum_RCM_approx)
!************************************************************************
!
! *sub_calcfg_l2_PO4* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Basis functions for solutes, case z > zbio
!
!               reac1, reac2        - mol./mol S released per organic carbon C
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: z, reac, Dtemp, ktemp, Qtemp , a2_M, alpha
REAL(8), INTENT(INOUT) :: e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, a2, b2
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(INOUT) :: Phi2
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: pfac
INTEGER :: i
REAL(8), DIMENSION(OS_RCM_fracs) :: loc_k

loc_k(:) = dum_RCM_approx(2,:)

a2 = (w-DSQRT(w**2.0D0+4.0D0*Dtemp*ktemp))/(2.0D0*Dtemp + OS_tolcte)
e = DEXP(z*a2)
dedz = a2*DEXP(z*a2)

b2=(w+DSQRT(w**2.0D0+4.0D0*Dtemp*ktemp))/(2.0D0*Dtemp + OS_tolcte)
f = DEXP(z*b2)
dfdz = b2*DEXP(z*b2)

pfac = 1.D0            !in fact, already has (1-por)/por

DO i=1, OS_RCM_fracs

   Phi2(i) = -pfac*loc_k(i)*(reac)*A22(i)/(Dtemp*aa2(i)**2.0D0-w*aa2(i)-ktemp + OS_tolcte)

END DO

g   =0.D0
dgdz=0.D0

DO i=1, OS_RCM_fracs
    g = g + Phi2(i)*DEXP(aa2(i)*z) 
    dgdz = Phi2(i)*aa2(i)*DEXP(aa2(i)*z)
END DO

IF (ktemp.NE.0) THEN        !CHECK: actually no need as ktemp always <> 0
   g =  g + Qtemp/(ktemp + OS_tolcte)   ! here problem if ktemp=0
END IF

IF (ktemp.NE.0) THEN            ! CHECK: think no need for this as always ktemp <> 0
   g = g + Qtemp/(ktemp + OS_tolcte)
END IF

IF (alpha.EQ.0) THEN            ! was z<=res.zox PO4 is independent of M (no info in alpha)
   p = 0.0D0
   dpdz = 0.0D0
   q = 0.0D0
   dqdz = 0.0D0
ELSE                        ! PO4 is dependent on M
   p = -alpha/(Dtemp*a2_M**2.0D0-w*a2_M-ktemp + OS_tolcte)*DEXP(a2_M*z)
   dpdz = a2_M*p
   q = 0.0D0
   dqdz = 0.0D0
END IF

END SUBROUTINE sub_calcfg_l2_PO4

!========================================================================

SUBROUTINE sub_calcfg_l1_M(z, Dtemp, ktemp, Qtemp, a1_P, b1_P, Phi1_P, alpha, e, dedz, &
                         & f, dfdz, g, dgdz, p, dpdz, q, dqdz, c1, d1)
!************************************************************************
!
! *sub_calcfg_l1_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Basis functions for solutes, case z <= zbio
!
!               reac1, reac2        - mol./mol S released per organic carbon C
!               depend1,   depend2 coming from other species
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!               and for dependent species
!               S(z) = A .* e(z) + B .* f(z) + C .* p(z) +  D.* q(z) + g(z)
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!    
REAL(8), INTENT(IN) :: z, Dtemp, ktemp, Qtemp, a1_P, b1_P, alpha                                
REAL(8), DIMENSION(3,OS_RCM_fracs), INTENT(INOUT) :: Phi1_P                                        
REAL(8), INTENT(INOUT) :: e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c1, d1
!
!*Local variables
!
INTEGER :: i

c1 = (w-DSQRT(w**2.0D0+4.0D0*Dtemp*ktemp))/(2.0D0*Dtemp + OS_tolcte)
p = DEXP(z*c1)
dpdz = c1*DEXP(z*c1)

d1 = (w+DSQRT(w**2.0D0+4.0D0*Dtemp*ktemp))/(2.0D0*Dtemp + OS_tolcte)
q = DEXP(z*d1)
dqdz = d1*DEXP(z*d1);

g    = 0.0D0
dgdz = 0.0D0

IF (alpha.NE.0.0D0) THEN      ! oxic layer: was z<=res.zox BUT problems at boundary. M is dependent on PO4
   c1 = 0.0D0
   d1 = 0.0D0
   e = -alpha/(Dtemp*a1_P**2.0D0-w*a1_P-ktemp + OS_tolcte)*DEXP(z*a1_P)
   dedz = a1_P*e
   f = -alpha/(Dtemp*b1_P**2.0D0-w*b1_P-ktemp + OS_tolcte)*DEXP(z*b1_P)
   dfdz = b1_P*f
   
   DO i=1, OS_RCM_fracs
      g = g + (Phi1_P(1,i)/(Dtemp*aa1(i)**2.0D0-w*aa1(i)-ktemp + OS_tolcte)*DEXP(z*aa1(i)) + &
             & Phi1_P(2,i)/(Dtemp*bb1(i)**2.0D0-w*bb1(i)-ktemp + OS_tolcte)*DEXP(z*bb1(i)) + &
             & Phi1_P(3,i)/(Dtemp*bb1(i)**2.0D0-w*bb1(i)-ktemp + OS_tolcte)*DEXP(z*bb1(i)))
      dgdz =  dgdz + (Phi1_P(1,i)/(Dtemp*aa1(i)**2.0D0-w*aa1(i)-ktemp + OS_tolcte)*DEXP(z*aa1(i))*aa1(i) + &
                    & Phi1_P(2,i)/(Dtemp*bb1(i)**2.0D0-w*bb1(i)-ktemp + OS_tolcte)*DEXP(z*bb1(i))*bb1(i) + &
                    & Phi1_P(3,i)/(Dtemp*bb1(i)**2.0D0-w*bb1(i)-ktemp + OS_tolcte)*DEXP(z*bb1(i))*bb1(i))
   END DO
   g = -alpha*g
   dgdz = -alpha*dgdz
    
ELSE                    ! anoxic layer: M is independent of PO4 (no value in alpha!)
   g = Qtemp/(ktemp + OS_tolcte)
   dgdz = 0.0D0
   e = 0.0D0
   dedz = 0.0D0
   f = 0.0D0
   dfdz = 0.0D0
END IF

END SUBROUTINE sub_calcfg_l1_M

!========================================================================

SUBROUTINE sub_calcfg_l2_M(z, ktemp, Qtemp, a2_P, b2_P, Phi2_P, alpha, e, dedz, f, dfdz, &
                         & g, dgdz, p, dpdz, q, dqdz, c2)
!************************************************************************
!
! *sub_calcfg_l1_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Basis functions for solutes, case z > zbio
!
!               reac1, reac2        - mol./mol S released per organic carbon C
!
!               General solution for solute S is given by
!               S(z) = A * e(z) + B * f(z) + g(z)
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!    
REAL(8), INTENT(IN) :: z, ktemp, Qtemp, a2_P, b2_P, alpha                                       
REAL(8), INTENT(INOUT) ::  e, dedz, f, dfdz, g, dgdz, p, dpdz, q, dqdz, c2            
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: Phi2_P                                       
!
!*Local variables
!
INTEGER :: i

c2 = 0

g    = 0.0D0
dgdz = 0.0D0

IF (alpha.NE.0.0D0) THEN            ! M is dependent of PO4, was z<=res.zox
   e = alpha/(w*a2_P + OS_tolcte)*DEXP(z*a2_P)
   dedz = a2_P*e
   f=alpha/(w*b2_P + OS_tolcte)*DEXP(z*b2_P)
   dfdz = b2_P*f
   p = 1.0D0                       ! DH CHECK/TODO: integration constant just C
   dpdz = 0.0D0
   q = 0.0D0
   dqdz = 0.0D0
   
   DO i=1, OS_RCM_fracs
     g = g + (Phi2_P(i)/(aa2(i) + OS_tolcte)*DEXP(aa2(i)*z))
     dgdz = dgdz +(Phi2_P(i)*DEXP(aa2(i)*z))
   END DO
   g = alpha/(w + OS_tolcte)*g
   dgdz = alpha/(w + OS_tolcte)*dgdz 
   
ELSE                        ! M is independent of PO4 - z > res.zox
   c2 = -ktemp/(w + OS_tolcte)
   p = DEXP(c2*z)
   dpdz = c2*DEXP(c2*z)
   q = 0.0D0
   dqdz = 0.0D0
   g = Qtemp/(ktemp + OS_tolcte)
   dgdz = 0.0D0
   e = 0.0D0
   dedz = 0.0D0
   f = 0.0D0
   dfdz = 0.0D0
END IF

END SUBROUTINE sub_calcfg_l2_M

!========================================================================

SUBROUTINE sub_prepfg_l12_PO4_M(reac1, ktempP, QtempP, zU, zL, D1P, D2P, alphaP, &
                              & ktempM, QtempM, D1M, D2M, alphaM, loc_mat_C, loc_vec_D, ltype, &
                              & dum_POC_conc_swi, dum_RCM_approx)
!************************************************************************
!
! *sub_prepfg_l12_PO4_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: reac1, ktempP, QtempP, zU, zL, D1P, D2P, alphaP, ktempM, QtempM, D1M, D2M, alphaM
REAL(8), DIMENSION(4,4), INTENT(INOUT) :: loc_mat_C
REAL(8), DIMENSION(1:4), INTENT(INOUT) :: loc_vec_D
INTEGER, INTENT(INOUT):: ltype
REAL(8), DIMENSION(OS_RCM_fracs), INTENT(IN) :: dum_POC_conc_swi
REAl(8), DIMENSION(2,OS_RCM_fracs), INTENT(IN) :: dum_RCM_approx
!
!*Local variables
!
REAL(8) :: e_zbio_l1_P, dedz_zbio_l1_P, f_zbio_l1_P, dfdz_zbio_l1_P, g_zbio_l1_P, dgdz_zbio_l1_P
REAL(8) :: p_zbio_l1_P, dpdz_zbio_l1_P, q_zbio_l1_P, dqdz_zbio_l1_P, a1_P, b1_P
REAL(8) :: e_zbio_l2_P, dedz_zbio_l2_P, f_zbio_l2_P, dfdz_zbio_l2_P, g_zbio_l2_P, dgdz_zbio_l2_P
REAL(8) :: p_zbio_l2_P, dpdz_zbio_l2_P, q_zbio_l2_P, dqdz_zbio_l2_P, a2_P, b2_P
REAL(8) :: e_zbio_l1_M, dedz_zbio_l1_M, f_zbio_l1_M, dfdz_zbio_l1_M, g_zbio_l1_M, dgdz_zbio_l1_M
REAL(8) :: p_zbio_l1_M, dpdz_zbio_l1_M, q_zbio_l1_M, dqdz_zbio_l1_M, a1_M, b1_M
REAL(8) :: e_zbio_l2_M, dedz_zbio_l2_M, f_zbio_l2_M, dfdz_zbio_l2_M, g_zbio_l2_M, dgdz_zbio_l2_M
REAL(8) :: p_zbio_l2_M, dpdz_zbio_l2_M, q_zbio_l2_M, dqdz_zbio_l2_M, a2_M
REAL(8) :: Vb, Fb
REAL(8), DIMENSION(3,OS_RCM_fracs) :: Phi1_P
REAL(8), DIMENSION(OS_RCM_fracs) :: Phi2_P
REAL(8), DIMENSION(1:4,1:4) :: mat_X, mat_Y
REAL(8), DIMENSION(1:4) ::  vec_Z
INTEGER :: loc_dim

Phi1_P = 0.D0 !(/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
Phi2_P = 0.D0 !(/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)

IF (zL.LE.zbio) THEN  ! wholly within bioturbated layer
   ltype = 1
ELSEIF (zU.GE.zbio) THEN ! wholly within non-bioturbated layer
   ltype = 2
ELSE             ! crossing boundary - sort out solution matching at zbio
   ltype = 3
   IF (zL.LE.zox) THEN       ! oxic layer -> call PO4 first
      CALL sub_calcfg_l1_PO4(z=zbio, reac=reac1, Dtemp=D1P, ktemp=ktempP, Qtemp=QtempP, &
                           & a1_M=0.0D0, b1_M=0.0D0, alpha=alphaP, e=e_zbio_l1_P, dedz=dedz_zbio_l1_P, f=f_zbio_l1_P, &
                           & dfdz=dfdz_zbio_l1_P, g=g_zbio_l1_P, dgdz=dgdz_zbio_l1_P, p=p_zbio_l1_P, dpdz=dpdz_zbio_l1_P, &
                           & q=q_zbio_l1_P, dqdz=dqdz_zbio_l1_P, a1=a1_P, b1=b1_P, Phi1=Phi1_P, &
                           & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
      CALL sub_calcfg_l2_PO4(z=zbio, reac=reac1, Dtemp=D2P, ktemp=ktempP, Qtemp=QtempP, &
                           & a2_M=0.0D0, alpha=alphaP, e=e_zbio_l2_P, dedz=dedz_zbio_l2_P, f=f_zbio_l2_P, &
                           & dfdz=dfdz_zbio_l2_P, g=g_zbio_l2_P, dgdz=dgdz_zbio_l2_P, p=p_zbio_l2_P, dpdz=dpdz_zbio_l2_P, &
                           & q=q_zbio_l2_P, dqdz=dqdz_zbio_l2_P, a2=a2_P, b2=b2_P, Phi2=Phi2_P, &
                           & dum_RCM_approx=dum_RCM_approx)
                           
      CALL sub_calcfg_l1_M(z=zbio, Dtemp=D1M, ktemp=ktempM, Qtemp=QtempM, a1_P=a1_P, b1_P=b1_P, &
                         & Phi1_P=Phi1_P, alpha=alphaM, e=e_zbio_l1_M, dedz=dedz_zbio_l1_M, &
                         & f=f_zbio_l1_M, dfdz=dfdz_zbio_l1_M, g=g_zbio_l1_M, dgdz=dgdz_zbio_l1_M, p=p_zbio_l1_M, &
                         & dpdz=dpdz_zbio_l1_M, q=q_zbio_l1_M, dqdz=dqdz_zbio_l1_M, c1=a1_M, d1=b1_M)
      CALL sub_calcfg_l2_M(z=zbio, ktemp=ktempM, Qtemp=QtempM, a2_P=a2_P, b2_P=b2_P, Phi2_P=Phi2_P, &
                         & alpha=alphaM, e=e_zbio_l2_M, dedz=dedz_zbio_l2_M, f=f_zbio_l2_M, dfdz=dfdz_zbio_l2_M, &
                         & g=g_zbio_l2_M, dgdz=dgdz_zbio_l2_M, p=p_zbio_l2_M, dpdz=dpdz_zbio_l2_M, q=q_zbio_l2_M, &
                         & dqdz=dqdz_zbio_l2_M, c2=a2_M)   

   ELSE              ! anoxic layer -> call M first
      CALL sub_calcfg_l1_M(z=zbio, Dtemp=D1M, ktemp=ktempM, Qtemp=QtempM, a1_P=0.0D0, b1_P=0.0D0, &
                         & Phi1_P=Phi1_P, alpha=alphaM, e=e_zbio_l1_M, dedz=dedz_zbio_l1_M, &
                         & f=f_zbio_l1_M, dfdz=dfdz_zbio_l1_M, g=g_zbio_l1_M, dgdz=dgdz_zbio_l1_M, p=p_zbio_l1_M, &
                         & dpdz=dpdz_zbio_l1_M, q=q_zbio_l1_M, dqdz=dqdz_zbio_l1_M, c1=a1_M, d1=b1_M)
      CALL sub_calcfg_l2_M(z=zbio, ktemp=ktempM, Qtemp=QtempM, a2_P=0.0D0, b2_P=0.0D0, Phi2_P=Phi2_P, &
                         & alpha=alphaM, e=e_zbio_l2_M, dedz=dedz_zbio_l2_M, f=f_zbio_l2_M, dfdz=dfdz_zbio_l2_M, &
                         & g=g_zbio_l2_M, dgdz=dgdz_zbio_l2_M, p=p_zbio_l2_M, dpdz=dpdz_zbio_l2_M, q=q_zbio_l2_M, &
                         & dqdz=dqdz_zbio_l2_M, c2=a2_M)                   
      CALL sub_calcfg_l1_PO4(z=zbio, reac=reac1, Dtemp=D1P, ktemp=ktempP, Qtemp=QtempP, &
                           & a1_M=a1_M, b1_M=b1_M, alpha=alphaP, e=e_zbio_l1_P, dedz=dedz_zbio_l1_P, f=f_zbio_l1_P, &
                           & dfdz=dfdz_zbio_l1_P, g=g_zbio_l1_P, dgdz=dgdz_zbio_l1_P, p=p_zbio_l1_P, dpdz=dpdz_zbio_l1_P, &
                           & q=q_zbio_l1_P, dqdz=dqdz_zbio_l1_P, a1=a1_P, b1=b1_P, Phi1=Phi1_P, &
                           & dum_POC_conc_swi=dum_POC_conc_swi, dum_RCM_approx=dum_RCM_approx)
      CALL sub_calcfg_l2_PO4(z=zbio, reac=reac1, Dtemp=D2P, ktemp=ktempP, Qtemp=QtempP, &
                           & a2_M=a2_M, alpha=alphaP, e=e_zbio_l2_P, dedz=dedz_zbio_l2_P, f=f_zbio_l2_P, &
                           & dfdz=dfdz_zbio_l2_P, g=g_zbio_l2_P, dgdz=dgdz_zbio_l2_P, p=p_zbio_l2_P, dpdz=dpdz_zbio_l2_P, &
                           & q=q_zbio_l2_P, dqdz=dqdz_zbio_l2_P, a2=a2_P, b2=b2_P, Phi2=Phi2_P, &
                           & dum_RCM_approx=dum_RCM_approx)
   END IF
   ! match solutions at zbio - continuous concentration and flux
   ! organize the data in matrices, and use the intrinsic fortran fct.
   ! DH: Maybe more efficient when written out !?
   !  |x1        |   | A_l |      | y1        | | A_r|    |z1|    always PO4 continuity
   !  |    .     |   | B_l |      |    .      | | B_r|    |z2|    always PO4 flux
   !  |      .   |   | C_l |   =  |      .    | | C_r|  + |z3|    always M continuity
   !  |       x16|   | D_l |      |        y16| | D_r|    |z4|    SD always M _diffusive_ flux  = 0 (cf org C)

   ! discontinuity constants
   Vb = 0.0D0
   Fb = 0.0D0

   ! weird FORTRAN matrices makes the transpose necessary
   ! matrix mat_X
   mat_X = TRANSPOSE(RESHAPE((/ e_zbio_l1_P, f_zbio_l1_P, p_zbio_l1_P, q_zbio_l1_P, &
                              & D1P*dedz_zbio_l1_P, D1P*dfdz_zbio_l1_P, D1P*dpdz_zbio_l1_P, &
                              & D1P*dqdz_zbio_l1_P, e_zbio_l1_M, f_zbio_l1_M, p_zbio_l1_M, &
                              & q_zbio_l1_M, D1M*dedz_zbio_l1_M, D1M*dfdz_zbio_l1_M, D1M*dpdz_zbio_l1_M, &
                              & D1M*dqdz_zbio_l1_M/), SHAPE(mat_X)))
   ! matrix mat_Y
   mat_Y = TRANSPOSE(RESHAPE((/ e_zbio_l2_P, f_zbio_l2_P, p_zbio_l2_P, q_zbio_l2_P, &
                              & D2P*dedz_zbio_l2_P, D2P*dfdz_zbio_l2_P, D2P*dpdz_zbio_l2_P, D2P*dqdz_zbio_l2_P, &
                              & e_zbio_l2_M, f_zbio_l2_M, p_zbio_l2_M, q_zbio_l2_M, D2M*dedz_zbio_l2_M, &
                              & D2M*dfdz_zbio_l2_M, D2M*dpdz_zbio_l2_M, D2M*dqdz_zbio_l2_M/), SHAPE(mat_Y)))
  
   vec_Z = (/ g_zbio_l2_P-g_zbio_l1_P + Vb, D2P*dgdz_zbio_l2_P - D1P*dgdz_zbio_l1_P + Fb - w*Vb, &
            & g_zbio_l2_M-g_zbio_l1_M + Vb, D2M*dgdz_zbio_l2_M - D1M*dgdz_zbio_l1_M + Fb - w*Vb /)

   loc_dim = 4
   CALL sub_matchsoln_PO4_M(mat_X, mat_Y, vec_Z, loc_dim, loc_mat_C, loc_vec_D)
END IF

END SUBROUTINE sub_prepfg_l12_PO4_M

!========================================================================

SUBROUTINE sub_matchsoln(E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l, &
                       & E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, &
                       & Vb, Db, ls_a, ls_b, ls_c, ls_d, ls_e, ls_f)
!************************************************************************
!
! *sub_matchsoln* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: E_l, F_l, G_l, dEdx_l, dFdx_l, dGdx_l
REAL(8), INTENT(IN) :: E_r, F_r, G_r, dEdx_r, dFdx_r, dGdx_r, Vb, Db
REAL(8), INTENT(INOUT) :: ls_a, ls_b, ls_c, ls_d, ls_e, ls_f
!
!*Local variables
!
REAL(8) :: alden, blden

! Match two solutions at a boundary:
! 'left' solution   y_l(x) = A_l*E_l(x) + B_l*F_l(x) + G_l(x)
! 'right' solution  y_r(x) = A_r*E_r(x) + B_r*F_l(x) + G_r(x)
!
! (Dis)continuity conditions at boundary:
!                   y_r(xb)    = y_l(xb)     + Vb
!                   dydx_r(xb) = dydx_l(xb)  + Db
!
! Find a,b,c,d,e,f such that:
!         | A_l |   =  | a  b | | A_r|  + |e|
!         | B_l |      | c  d | | B_r|    |f|

alden = dFdx_l*E_l - F_l*dEdx_l + OS_tolcte
ls_a     = (dFdx_l*E_r - F_l*dEdx_r)/alden
ls_b     = (dFdx_l*F_r - F_l*dFdx_r)/alden
ls_e     = (F_l*(dGdx_l - dGdx_r + Db) + dFdx_l*(-G_l + G_r - Vb))/alden

blden = dEdx_l*F_l - E_l*dFdx_l + OS_tolcte
ls_c     = (dEdx_l*E_r - E_l*dEdx_r)/blden
ls_d     = (dEdx_l*F_r - E_l*dFdx_r)/blden;
ls_f     = (E_l*(dGdx_l - dGdx_r + Db) + dEdx_l*(-G_l+G_r - Vb))/blden

END SUBROUTINE sub_matchsoln

!========================================================================

SUBROUTINE sub_xformsoln(E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f, Etr, Ftr, Gtr, dEtdx, dFtdx, dGtdx)
!************************************************************************
!
! *sub_xformsoln* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: E, F, G, dEdx, dFdx, dGdx, ls_a , ls_b , ls_c , ls_d , ls_e ,ls_f
REAL(8), INTENT(INOUT) :: Etr, Ftr, Gtr, dEtdx, dFtdx, dGtdx
!
!*Local variables
!

! Find 'transformed' soln such that in layer l,
!    y_l = A_r*et + B_r*ft + gt
! (ie l soln written in terms of r solution coefficents A_r, B_r)

Etr   = ls_a*E    + ls_c*F
dEtdx = ls_a*dEdx + ls_c*dFdx
Ftr   = ls_b*E    + ls_d*F
dFtdx = ls_b*dEdx + ls_d*dFdx
Gtr   = G         + ls_e*E      + ls_f*F
dGtdx = dGdx      + ls_e*dEdx   + ls_f*dFdx

END SUBROUTINE sub_xformsoln

!========================================================================

SUBROUTINE sub_solve2eqn(a, b, c, d, e, f, x, y)
!************************************************************************
!
! *sub_solve2eqn* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -  Find soln of
!                | a    b |  |x|   = | e |
!                | c    d |  |y|     | f |
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), INTENT(IN) :: a, b, c, d, e, f
REAL(8), INTENT(OUT) :: x, y
!
!*Local variables
!
REAL(8) :: det

det = a*d-b*c+OS_tolcte
IF (det.EQ.0.0D0) THEN
    print*,'det too small ', det
END IF
x = (e*d-b*f)/det
y = (a*f-e*c)/det

END SUBROUTINE sub_solve2eqn

!========================================================================

SUBROUTINE sub_matchsoln_PO4_M(dum_mat_X, dum_mat_Y, dum_vec_Z, dum_dim, loc_mat_C, loc_vec_D)
!************************************************************************
!
! *sub_matchsoln_PO4_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: dum_dim                      ! dimension of the matrices
REAL(8), DIMENSION(1:dum_dim,1:dum_dim), INTENT(IN) :: dum_mat_X, dum_mat_Y
REAL(8), DIMENSION(1:dum_dim), INTENT(IN) ::  dum_vec_Z
REAL(8), DIMENSION(1:dum_dim,1:dum_dim), INTENT(INOUT) :: loc_mat_C
REAL(8), DIMENSION(1:dum_dim), INTENT(INOUT) ::  loc_vec_D
!
!*Local variables
!
REAL(8), DIMENSION(1:dum_dim,1:dum_dim) :: loc_mat_X, loc_inv_mat_X

! Match four solutions at a boundary:
!  for PO4
! 'left' solution   y_l(z) = A_l*E_l(z) + B_l*F_l(z) + C_l*P_l(z) + D_l*Q_l(z) + G_l(z)
! 'right' solution  y_r(z) = A_r*E_r(z) + B_r*F_r(z) + C_r*P_r(z) + D_r*Q_r(z) + G_r(z)
!
!  and the same for M
!
! (Dis)continuity conditions at boundary:
!                   y_r(xb)    = y_l(xb)     + Vb
!                   dydx_r(xb) = dydx_l(xb)  + Db
!         | A_l |         | A_r|
!         | B_l |         | B_r|
!     X   | C_l |   =  Y  | C_r|  + Z
!         | D_l |         | D_r|
!
! Find C and D such that:
!         | A_l |         | A_r|
!         | B_l |         | B_r|
!         | C_l |   =  C  | C_r|  +  D
!         | D_l |         | D_r|

! save matrix locally, as the original matrix loc_mat_X(4,4) will be destroyed during the calculation
loc_mat_X = dum_mat_X
! calculate loc_mat_X^{-1}
CALL sub_inverse(loc_mat_X,loc_inv_mat_X,dum_dim)
loc_mat_C = MATMUL(loc_inv_mat_X, dum_mat_Y)
loc_vec_D = MATMUL(loc_inv_mat_X, dum_vec_Z)

END SUBROUTINE sub_matchsoln_PO4_M

!========================================================================

SUBROUTINE sub_xformsoln_PO4_M(dum_EFPQ_P, dum_EFPQ_M, dum_dEFPQdz_P, dum_dEFPQdz_M, &
                             & dum_g_P, dum_g_M, dum_dgdz_P, dum_dgdz_M, dum_mat_C, dum_vec_D, &
                             & loc_EFPQ_P_t, loc_G_P_t, loc_dEFPQ_P_t, loc_dG_P_t, loc_EFPQ_M_t, loc_G_M_t, &
                             & loc_dEFPQ_M_t, loc_dG_M_t)
!************************************************************************
!
! *sub_xformsoln_PO4_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - 
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8), DIMENSION(1:4), INTENT(IN) :: dum_EFPQ_P, dum_dEFPQdz_P, dum_EFPQ_M, dum_dEFPQdz_M
REAL(8), INTENT(IN) :: dum_g_P, dum_g_M, dum_dgdz_P, dum_dgdz_M
REAL(8), DIMENSION(4, 4), INTENT(IN) :: dum_mat_C
REAL(8), DIMENSION(1:4), INTENT(IN) ::  dum_vec_D
REAL(8), DIMENSION(1:4), INTENT(INOUT) :: loc_EFPQ_P_t, loc_dEFPQ_P_t, loc_EFPQ_M_t, loc_dEFPQ_M_t
REAL(8), INTENT(INOUT) :: loc_G_P_t, loc_dG_P_t, loc_G_M_t, loc_dG_M_t
!
!*Local variables
!

! Find 'transformed' soln such that in layer l,
!    y_l = A_r*et + B_r*ft + gt
! (ie l soln written in terms of r solution coefficents A_r, B_r)
!
! here save values in matrices - as this saves a lot of code

loc_EFPQ_P_t = MATMUL(TRANSPOSE(dum_mat_C), dum_EFPQ_P)     ! DH TODO: check multiplication with transpose, especially if vector*vector
loc_dEFPQ_P_t = MATMUL(TRANSPOSE(dum_mat_C), dum_dEFPQdz_P)
loc_G_P_t = DOT_PRODUCT(dum_vec_D, dum_EFPQ_P)+dum_g_P
loc_dG_P_t = DOT_PRODUCT(dum_vec_D,dum_dEFPQdz_P)+dum_dgdz_P

loc_EFPQ_M_t = MATMUL(TRANSPOSE(dum_mat_C), dum_EFPQ_M)
loc_dEFPQ_M_t = MATMUL(TRANSPOSE(dum_mat_C), dum_dEFPQdz_M)
loc_G_M_t = DOT_PRODUCT(dum_vec_D, dum_EFPQ_M) + dum_g_M
loc_dG_M_t = DOT_PRODUCT(dum_vec_D,dum_dEFPQdz_M) + dum_dgdz_M

END SUBROUTINE sub_xformsoln_PO4_M

!========================================================================

SUBROUTINE sub_solve2eqn_PO4_M(dum_mat_X, dum_vec_Y, dum_A, dum_B, dum_C, dum_dim)
!************************************************************************
!
! *sub_xformsoln_PO4_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description -  Find solution of
!                | x1  .  . x4|  |A|     | y1 |
!                |     .      |  |B|     | y2 |
!                |       .    |  |C|   = | y3 |
!                | .       x16|  |D|     | y4 |
!
! Reference -
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: dum_dim                             ! dim of input matrix to inverse
REAL(8), DIMENSION(1:dum_dim,1:dum_dim), INTENT(IN) :: dum_mat_X
REAL(8), DIMENSION(1:dum_dim), INTENT(IN) :: dum_vec_Y
REAL(8), INTENT(INOUT) :: dum_A, dum_B, dum_C
!
!*Local variables
!
REAL(8), DIMENSION(1:dum_dim,1:dum_dim) :: loc_mat_X, loc_inv_mat_X
REAL(8), DIMENSION(1:dum_dim) :: loc_vec_Z

! save matrix locally, as the original matrix dum_mat_X(4,4) will be destroyed during the calculation
loc_mat_X = dum_mat_X
! calculate loc_mat_X^{-1}
CALL sub_inverse(loc_mat_X,loc_inv_mat_X,dum_dim)
loc_vec_Z = MATMUL(loc_inv_mat_X, dum_vec_Y)

dum_A = loc_vec_Z(1)
dum_B = loc_vec_Z(2)
dum_C = loc_vec_Z(3)

END SUBROUTINE sub_solve2eqn_PO4_M

!========================================================================

SUBROUTINE sub_inverse(a, c,n)
!************************************************************************
!
! *sub_xformsoln_PO4_M* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - Inverse matrix
!               input ...
!                     a(n,n) - array of coefficients for matrix A
!                     n      - dimension
!               output ...
!                     c(n,n) - inverse matrix of A
!               comments ...
!                     the original matrix a(n,n) will be destroyed
!                     during the calculation
!
! Reference - Based on Doolittle LU factorization for Ax=b
!             Alex G. December 2009
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: n
REAL(8), DIMENSION(1:n,1:n), INTENT(INOUT) :: a, c
!
!*Local variables
!
REAL(8), DIMENSION(1:n,1:n) :: L, U
REAL(8), DIMENSION(1:n) :: b, d, x
REAL(8) :: coeff
INTEGER :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 allows such operations on matrices
L = 0.0D0
U = 0.0D0
b = 0.0D0

! step 1: forward elimination
DO k=1, n-1
  DO i=k+1,n
    coeff=a(i,k)/(a(k,k) + OS_tolcte)
    L(i,k) = coeff
    DO j=k+1,n
      a(i,j) = a(i,j)-coeff*a(k,j)
    END DO
  END DO
END DO

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
DO i=1,n
  L(i,i) = 1.0D0
END DO
! U matrix is the upper triangular part of A
DO j=1,n
  DO i=1,j
    U(i,j) = a(i,j)
  END DO
END DO

! Step 3: compute columns of the inverse matrix C
DO k=1,n
  b(k)=1.0D0
  d(1) = b(1)
  ! Step 3a: Solve Ld=b using the forward substitution
  DO i=2,n
    d(i)=b(i)
    DO j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    END DO
  END DO
  ! Step 3b: Solve Ux=d using the back substitution
  ! DH: check for division by zero
  x(n)=d(n)/(U(n,n)+OS_tolcte)
  DO i = n-1,1,-1
    x(i) = d(i)
    DO j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    END DO
    x(i) = x(i)/(u(i,i) + OS_tolcte)
  END DO
  ! Step 3c: fill the solutions x(n) into column k of C
  DO i=1,n
    c(i,k) = x(i)
  END DO
  b(k)=0.0D0
END DO

END SUBROUTINE sub_inverse

!========================================================================

FUNCTION FUN_zbrent(func,x1,x2,tol)
!************************************************************************
!
! *FUN_zbrent* calculate  
!
! Author - Dominik Hülse and Sebastiaan van de Velde
!
! Version - @(OMENSED)OMENSED_module.f90  V0.1
!
! Description - calculate root of func in the interval [x1,x2]
!
! Reference - 
!
! Calling program - 
!
!************************************************************************    
         
IMPLICIT NONE
!
!*Arguments
!
REAL(8) :: FUN_zbrent,tol,x1,x2,func
!
!*Local variables
!
INTEGER :: ITMAX
REAL(8) :: EPS
EXTERNAL func
PARAMETER (ITMAX=100,EPS=3.e-8)
INTEGER :: iter
REAL(8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        
a=x1
b=x2
fa=func(a)
fb=func(b)
IF ((fa.GT.0..AND.fb.GT.0.).OR.(fa.LT.0..AND.fb.LT.0.)) THEN
   print*,'root must be bracketed for FUN_zbrent'
   STOP
ELSE
   c=b
   fc=fb
   DO 11 iter=1,ITMAX
     IF ((fb.GT.0..AND.fc.GT.0.).OR.(fb.LT.0..AND.fc.LT.0.)) THEN
        c=a
        fc=fa
        d=b-a
        e=d
     END IF
     IF (DABS(fc).LT.DABS(fb)) THEN
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     END IF
        tol1=2.D0*EPS*DABS(b)+0.5D0*tol
        xm=.5*(c-b)
     IF (DABS(xm).LE.tol1.OR.fb.EQ.0.) THEN
         FUN_zbrent=b
         RETURN
     END IF
     IF (DABS(e).GE.tol1.AND.DABS(fa).GT.DABS(fb)) THEN
        s=fb/fa
        IF (a.EQ.c) THEN
           p=2.D0*xm*s
           q=1.D0-s
        ELSE
           q=fa/fc
           r=fb/fc
           p=s*(2.D0*xm*q*(q-r)-(b-a)*(r-1.D0))
           q=(q-1.D0)*(r-1.D0)*(s-1.D0)
        END IF
        IF (p.GT.0.) q=-q
        p=DABS(p)
        IF(2.D0*p.LT.MIN(3.D0*xm*q-DABS(tol1*q),DABS(e*q))) THEN
           e=d
           d=p/q
        ELSE
           d=xm
           e=d
        END IF
      ELSE
        d=xm
        e=d
      END IF
      a=b
      fa=fb
      IF (DABS(d) .GT. tol1) THEN
         b=b+d
      ELSE
         b=b+SIGN(tol1,xm)
      END IF
      fb=func(b)
11    CONTINUE
END IF
print*,'FUN_zbrent exceeding maximum iterations'
FUN_zbrent=b
STOP

RETURN

END FUNCTION FUN_zbrent

!========================================================================


END MODULE OMEN_Functions
