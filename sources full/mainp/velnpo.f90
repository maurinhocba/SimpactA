 SUBROUTINE velnpo (npoin, nvelr, ifpre, nesdf, npsdf,      &
                    ftsdf, velor, veloc, velnp)
 !***********************************************************************
 !
 !     rearrange velocity vector
 !
 !***********************************************************************
 !$ USE omp_lib
 USE ctrl_db, ONLY : nrotd,vefac
 USE kinc_db, ONLY : nn,nm
 IMPLICIT NONE

 INTEGER (kind=4), INTENT(IN) :: npoin, nvelr, ifpre(:,:),  &
                                 nesdf(:),npsdf(:)
 REAL (kind=8), INTENT(IN)  :: ftsdf(:), veloc(:), velor(:,:)
 REAL (kind=8), INTENT(OUT) :: velnp(:,:)

 ! local
 INTEGER (kind=4) :: ieq,ipoin,idofn,iq,ib,ie,i,nv1
 REAL (kind=8) :: auxil,vf1

 vf1 = vefac + 1d0
 nv1 = nvelr+1
 !$OMP DO  &
 !$OMP PRIVATE(ipoin,idofn,ieq,ib,ie,auxil,i,iq)
 DO ipoin = 1,npoin
   DO idofn = 1,nrotd
     ieq = ifpre(idofn,ipoin)
     SELECT CASE (ieq)
     CASE (1:)                                 !active dof
       velnp(idofn,ipoin) = (veloc(ieq) + velnp(idofn,ipoin) * vefac )/vf1
     CASE (0)                                  !null dof
       velnp(idofn,ipoin) = 0.
     CASE ( -nm:-1)                            !slave DOFs
       ib = npsdf(-ieq)
       ie = npsdf(-ieq+1)-1
       auxil = 0d0
       DO i = ib,ie
         iq = nesdf(i)
         IF(iq > 0)THEN
           auxil = auxil + veloc(iq)*ftsdf(i)
         ELSE IF(iq < -nn)THEN
           auxil = auxil + velor(-iq-nn,nv1)*ftsdf(i)
         END IF
       END DO
       velnp(idofn,ipoin) = ( auxil + velnp(idofn,ipoin) * vefac )/vf1
     CASE ( :-nn)                           !fixed DOFs
       velnp(idofn,ipoin) = velor(-ieq-nn,nv1)
     END SELECT
   END DO
 END DO
 !$OMP END DO

 RETURN
 END SUBROUTINE velnpo
