 SUBROUTINE coract(lambd,fixvl,ddisp,coora,euler)

 !     UPDATES CONFIGURATION AND NODAL SYSTEMS
 !$ USE omp_lib
 USE ctrl_db, ONLY : ndime,npoin,ndofn,neulr,nrotd
 USE kinc_db, ONLY : nn,nn1,nm,nescv,ndepd,naris,&
                     nndpd,distd,npsdf,nesdf,ftsdf
 USE npo_db, ONLY : naeul,ifpre,psia
 USE rfdp_db, ONLY : nrfdp
 USE nsld_db, ONLY : nnsld
 IMPLICIT NONE
 REAL (kind=8),INTENT(IN) :: ddisp(:),lambd,fixvl(:)
 REAL (kind=8),INTENT(IN OUT) :: coora(:,:)
 REAL (kind=8), POINTER :: euler(:,:)

 INTEGER (kind=4) :: ipoin,idofn,ieq,i,ie,ib,iq
 REAL    (kind=8) :: theta(3),auxil
 !REAL(kind=8) :: ls(9,4)
 ! REAL (kind=8),PARAMETER :: pi = 3.1415926535898_8

 INTERFACE
   INCLUDE 'actrot.h'
   INCLUDE 'arisfi.h'
   INCLUDE 'fixdpd.h'
   INCLUDE 'fixrfd.h'
   INCLUDE 'fixsld.h'
 END INTERFACE

 !     translational degrees of freedom
 !$OMP DO  PRIVATE(ipoin,idofn,ieq,ib,ie,auxil,i,iq)
 DO ipoin=1,npoin
   DO idofn=1,ndime
     ieq = ifpre(idofn,ipoin)
     SELECT CASE (ieq)
     CASE (:-nn1)          !fixed DOF
       coora(idofn,ipoin) = coora(idofn,ipoin)+fixvl(-ieq-nn)*lambd
     CASE (-nn)
       !nothing
     CASE (-nm:-1)          !slave DOF
       ib = npsdf(-ieq)
       ie = npsdf(-ieq+1)-1
       auxil = 0d0
       DO i = ib,ie
         iq = nesdf(i)
         IF(iq > 0)THEN
           auxil = auxil + ddisp(iq)*ftsdf(i)
         ELSE IF(iq < -nn)THEN
           auxil = auxil + fixvl(-iq-nn)*ftsdf(i)*lambd
         END IF
       END DO
       coora(idofn,ipoin) = coora(idofn,ipoin) + auxil
     CASE (0)                 !null DOF
       !nothing
     CASE (1:)                !active DOF
       coora(idofn,ipoin) = coora(idofn,ipoin) + ddisp(ieq)
     END SELECT
   END DO
 END DO
 !$OMP END DO
 ! rotational degrees of freedom
 IF (neulr == 1) THEN
   idofn=ndime+1
   !$OMP DO  &
   !$OMP PRIVATE(ipoin,ieq,ib,ie,auxil,i,iq)
   DO ipoin=1,npoin
     IF( .NOT. naeul(ipoin) )CYCLE
     ieq = ifpre(idofn,ipoin)
     SELECT CASE (ieq)
     CASE (:-nn1)                !fixed DOF
       euler(1,ipoin) = euler(1,ipoin)+lambd*fixvl(-ieq-nn)
     CASE (-nn:-1)               !slave DOF
       ib = npsdf(-ieq)
       ie = npsdf(-ieq+1)-1
       auxil = 0d0
       DO i = ib,ie
         iq = nesdf(i)
         IF(iq > 0)THEN
           auxil = auxil + ddisp(iq)*ftsdf(i)
         ELSE IF(iq < -nn)THEN
           auxil = auxil + fixvl(-iq-nn)*ftsdf(i)*lambd
         END IF
       END DO
       euler(1,ipoin) = euler(1,ipoin)+auxil
     !CASE (0)                     !null DOF
     !  euler(1,ipoin) = euler(1,ipoin)
     CASE (1:)                     !active  DOF
       euler(1,ipoin) = euler(1,ipoin)+ddisp(ieq)
     END SELECT
   END DO
   !$OMP END DO
 ELSE IF (neulr == 9) THEN
   theta(3) = 0d0
   !ls(:,:) = euler(:,(/1,3,23,43/))
   !$OMP DO  &
   !$OMP PRIVATE(ipoin,idofn,ieq,ib,ie,auxil,i,iq,theta)
   DO ipoin=1,npoin
     IF( .NOT. naeul(ipoin) )CYCLE
     DO idofn=ndime+1,nrotd
       ieq = ifpre(idofn,ipoin)
       SELECT CASE (ieq)
       CASE (:-nn1)        !fixed DOF
         theta(idofn-ndime) = lambd*fixvl(-ieq-nn)
       CASE (-nn:-1)        !slave DOF
         ib = npsdf(-ieq)
         ie = npsdf(-ieq+1)-1
         auxil = 0d0
         DO i = ib,ie
           iq = nesdf(i)
           IF(iq > 0)THEN
             auxil = auxil + ddisp(iq)*ftsdf(i)
           ELSE IF(iq < -nn)THEN
            auxil = auxil + fixvl(-iq-nn)*ftsdf(i)*lambd
           END IF
         END DO
         theta(idofn-ndime) = auxil
       CASE (0)               !null DOF
         theta(idofn-ndime) = 0d0
       CASE (1:)              !active or DOF
         theta(idofn-ndime) = ddisp(ieq)
       END SELECT
     END DO
     CALL actrot(euler(1:9,ipoin),theta)
   END DO
   !$OMP END DO
 END IF
 IF(ndofn == 8) THEN
   DO ipoin=1,npoin
     DO idofn=1,2
       ieq = ifpre(idofn+6,ipoin)
       SELECT CASE (ieq)
       CASE (:-nn1)        !fixed DOF
         psia(idofn,ipoin) = lambd*fixvl(-ieq-nn)
       CASE (-nn:-1)        !slave DOF
         ib = npsdf(-ieq)
         ie = npsdf(-ieq+1)-1
         auxil = 0d0
         DO i = ib,ie
           iq = nesdf(i)
           IF(iq > 0)THEN
             auxil = auxil + ddisp(iq)*ftsdf(i)
           ELSE IF(iq < -nn)THEN
            auxil = auxil + fixvl(-iq-nn)*ftsdf(i)*lambd
           END IF
         END DO
         psia(idofn,ipoin) = auxil
       CASE (0)               !null DOF
         !nothing
       CASE (1:)              !active or DOF
         psia(idofn,ipoin) = psia(idofn,ipoin) + ddisp(ieq)
       END SELECT
     END DO
   END DO
 END IF

 IF(ndepd > 0) CALL fixdpd(coora,euler)
 !WRITE(58,"(3(3e20.12,/))")euler(:,(/1,3,23,43/))-ls
 IF(nrfdp > 0) CALL fixrfd(coora)
 IF(naris > 0) CALL arisfi(naris,coora,nndpd(1:3,ndepd+1:ndepd+naris),euler)
 IF(nnsld > 0) CALL fixsld(coora,euler)

 RETURN

 END SUBROUTINE coract
