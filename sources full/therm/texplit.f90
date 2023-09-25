 SUBROUTINE texplit( )
 !********************************************************************
 !
 ! *** time stepping routine
 !
 !********************************************************************
 USE ctrl_db, ONLY: tdtime, ttime, ndoft, nheat, npoin, therm, neqt
 USE curv_db
 !      USE kinc_db

 USE ift_db, ONLY : npret,lctmp,nprev,prtmp,fxtem
 USE npo_db, ONLY : tempe,dtemp,tresi,heata,heati,iftmp,trate,tmass
 USE heat_db, ONLY : ncsur
 IMPLICIT NONE
! local variables
 INTEGER (kind=4), PARAMETER :: nn = 1000000, nn1 = 1000001
 INTEGER (kind=4) :: ieq,i,n1,nt0,j
 REAL (kind=8) :: dt2,facts(nheat),factv(npret)
 REAL (kind=8), SAVE :: tdtold
! Functions
 REAL (kind=8):: functs

 !   compute prescribed temperatures
 DO i=1,npret
   factv(i) = functs (lctmp(i),ttime)*prtmp(nprev+1,i)
 END DO

 IF(npret > 0)THEN
   DO i=1,nprev
     prtmp(i,npret+1)= DOT_PRODUCT(prtmp(i,1:npret),factv(1:npret))
   END DO
 END IF

 dtemp = tempe                      !keep original values

 IF( therm )THEN                    !for heat transfer problems

   IF(ttime == 0d0) tdtold = 0d0    !first time
   dt2 = (tdtime+tdtold)/2d0        !average time increment

   IF( ncsur > 0 ) CALL tcsurf( )    !compute convection and radiation terms

   trate = -tresi                   !initializes residual

   !   compute and add external equivalent heats
   DO i=1,nheat
     facts(i) = functs (heata(i),ttime)*heati(neqt+1,i)
   END DO
   IF(nheat > 0)THEN
     n1  = nheat+1
     DO ieq=1,neqt
       heati(ieq,n1) = DOT_PRODUCT(heati(ieq,1:nheat),facts(1:nheat))
       trate(ieq) = trate(ieq) + heati(ieq,n1)
     END DO
   END IF
   !   integrate equations
   DO ieq = 1,neqt
     trate(ieq) =  trate(ieq)/tmass(ieq)*dt2      !increment in fact
   END DO

   !  update temperatures

   DO i=1,npoin             !loop over each node
     DO j=1,ndoft
       nt0 = iftmp(j,i)       !DOF
       SELECT CASE (nt0)
       CASE (1: )             ! active DOF
         tempe(j,i) = tempe(j,i) + trate(nt0)
       CASE (0 )
         !nothing
       CASE (-nn:-1)          ! fixed DOF
         tempe(j,i) = fxtem(-nt0)
       CASE (:-nn1)           ! prescribed DOF
         tempe(j,i) = prtmp(-nt0-nn,npret+1)
       END SELECT
     END DO
   END DO

   tdtold = tdtime

 ELSE  ! itemp only

   !  update temperatures

   DO i=1,npoin             !loop over each node
     DO j=1,ndoft             !loop over each node
       nt0 = iftmp(j,i)       !DOF
       SELECT CASE (nt0)
       CASE (0 )
         !nothing
       CASE (-nn:-1)
         tempe(j,i) = fxtem(-nt0)
       CASE (:-nn1)
         tempe(j,i) = prtmp(-nt0-nn,npret+1)
       END SELECT
     END DO
   END DO

 END IF

 dtemp = tempe - dtemp  !incremental temperature

 RETURN
 END SUBROUTINE texplit
