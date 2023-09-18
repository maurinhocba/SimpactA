 SUBROUTINE rigbdc ( nndp )
 !
 !     This routine computes:
 !          nodal dependencies of nodes in rigid bodies
 !
 USE ctrl_db, ONLY : npoin,ndime,neulr
 USE kinc_db, ONLY : nndpd,distd,npsdf,nesdf,ftsdf,ndumn
 USE npo_db,  ONLY : coora,ifpre,euler,naeul
 IMPLICIT NONE
 INTEGER (kind=4), INTENT(IN) :: nndp  ! other dependent nodes

 INTEGER (kind=4) i,j,k,l,m,mm,n,n1,n2
 REAL (kind=8) :: d(ndime),t1,t2
 LOGICAL eul

 eul = neulr > 0

 DO n=nndp+1,nndp+ndumn
   n1 = nndpd(1,n)  !slave node label
   n2 = nndpd(2,n)  !master node label
   d = coora(1:ndime,n1) - coora(1:ndime,n2)    !distance vector between nodes (master->slave)
   IF( eul )THEN
     naeul(n2) = .TRUE.             !local system required at master node
     IF(ndime == 3) THEN            !for 3-D problems
       distd(1:3,n) = MATMUL(d,RESHAPE(euler(1:9,n2),(/3,3/))) !material components
       ! translational degrees of freedom
       j = - ifpre(1,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k+1) = -euler(4,n2)*distd(3,n) + euler(7,n2)*distd(2,n)
       ftsdf(k+2) = -euler(7,n2)*distd(1,n) + euler(1,n2)*distd(3,n)
       ftsdf(k+3) = -euler(1,n2)*distd(2,n) + euler(4,n2)*distd(1,n)

       j = - ifpre(2,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k+1) = -euler(5,n2)*distd(3,n) + euler(8,n2)*distd(2,n)
       ftsdf(k+2) = -euler(8,n2)*distd(1,n) + euler(2,n2)*distd(3,n)
       ftsdf(k+3) = -euler(2,n2)*distd(2,n) + euler(5,n2)*distd(1,n)

       j = - ifpre(3,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k+1) = -euler(6,n2)*distd(3,n) + euler(9,n2)*distd(2,n)
       ftsdf(k+2) = -euler(9,n2)*distd(1,n) + euler(3,n2)*distd(3,n)
       ftsdf(k+3) = -euler(3,n2)*distd(2,n) + euler(6,n2)*distd(1,n)
       ! IF slave node have local system
       IF( naeul(n1) )THEN   ! rotational degrees of freedom
         m=1                        !first position in rotation matrix
         DO i=1,ndime                 !for each rotational DOF
           j = - ifpre(3+i,n1)        !slave DOF order
           k = npsdf(j)               !first position for this DOF in NESDF & FTSDF
           mm = 1
           DO l=1,ndime
             ftsdf(k) = DOT_PRODUCT(euler(m:m+2,n1),euler(mm:mm+2,n2)) !factor
             k = k+1                  !counter in NESDF & FTSDF
             mm = mm + 3
           END DO
           m = m+3
         END DO
       END IF
     ELSE ! ndime = 2 (2-D problems)
       t1 = COS(euler(1,n2))                !local system of master node
       t2 = SIN(euler(1,n2))
       distd(1,n)= t1*d(1)+t2*d(2)          !material components
       distd(2,n)=-t2*d(1)+t1*d(2)          !of distance vector

       j = - ifpre(1,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k+1) = -t2*distd(1,n) - t1*distd(2,n) !factor for rotational DOF

       j = - ifpre(2,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k+1) = t1*distd(1,n) - t2*distd(2,n) !factor for rotational DOF

     END IF
   ELSE  !model does not include rotational DOFs
     distd(1:ndime,n) = d(1:ndime) !material components
   END IF
 END DO
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rigbdc
