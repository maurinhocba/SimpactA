 SUBROUTINE readpd(npsdf,nesdf,ftsdf,j,np)
 !**********************************************************************
 !
 !     APPLIES and generates DATA for dependant nodes
 !
 !**********************************************************************
 USE lispa0, ONLY : lures
 USE ndp_db
 USE ctrl_db, ONLY : ndime,ndofn,nrotd,npoin,neulr
 USE npo_db,  ONLY : naeul,ifpre,coora,euler,label
 USE kinc_db, ONLY : ndepd,nn,distd,nndpd
 USE outp_db, ONLY: iwrit

 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: np
 INTEGER (kind=4),INTENT(IN OUT) :: j,npsdf(np),nesdf(4*np)
 REAL (kind=8), INTENT(IN OUT) :: ftsdf(4*np)

 INTEGER (kind=4) :: i,k,l,m,n,n1,n2,ll
 REAL (kind=8) :: d(ndime),t1,t2,ang(3)
 LOGICAL :: eul

 INTEGER(kind=4) chnode
 INTERFACE
   INCLUDE 'inrotm.h'
  END INTERFACE

 eul = neulr > 0
 IF((ndime == 3 .AND. ndofn < 6).OR.(ndime == 2 .AND. ndofn < 3).AND. iwrit == 1)&
    WRITE(lures,"(' its not posible to have true dependent nodes in',i2,&
          &    '-D with only',i2,' DOF')",ERR=9999) ndime,ndofn

 IF (iwrit == 1 .AND. nndp > 0 ) WRITE(lures,"(//,3x,'Slave/Master nodes Data',//&
                 & '   Slave Node     Master Node')",ERR=9999)

 DO n=1,nndp          !loop for each dependant node

   n1 =ndpdat(1,n)    !slave node label
   n2 =ndpdat(2,n)    !master node label
   IF (iwrit == 1) WRITE(lures, '(2i10)',ERR=9999) n1, n2   !ECHO
   n1 = chnode(n1)              !slave internal node
   n2 = chnode(n2)              !master internal node
   nndpd(1,n) = n1              !copy slave into array
   nndpd(2,n) = n2              !copy master into array
   DO i = 1,nrotd               !for each (slave) DOF
     IF(ifpre(i,n1) < 0 ) THEN  !check that the DOF is not already slave
       WRITE(lures,"(' DOF already declared Slave:  ',&
                  &  'Node',i6,' DOF',i2)",ERR=9999) label(n1) ,i
       CALL runend('READPD: Inconsistent Input Data    ')
     END IF
   END DO
   d = coora(1:ndime,n1) - coora(1:ndime,n2)    !distance vector between nodes (master->slave)
   IF( eul )THEN
     IF( .NOT. naeul(n2) )THEN
       naeul(n2) = .TRUE.             !local system required at master node
       IF(ndime == 3) THEN            !for 3-D problems
         IF( ALL(euler(4:9,n2) == 0d0))THEN
           ang = euler(1:3,n2)              !Recover the Euler angles (radians)
           CALL inrotm(ang,euler(1:9,n2))   !compute rotation matrix
         END IF
       END IF
     END IF
     IF(ndime == 3) THEN            !for 3-D problems
       distd(1:3,n) = MATMUL(d,RESHAPE(euler(1:9,n2),(/3,3/))) !material components
       ! translational degrees of freedom
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(1,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1 + ndime   !compute next pointer
       nesdf(k)   = 1+10*n2         !store master node and dof until later
       nesdf(k+1) = 4+10*n2         !
       nesdf(k+2) = 5+10*n2         !
       nesdf(k+3) = 6+10*n2         !
       ftsdf(k)   = 1d0             !factor for translational DOF & rotational DOFs
       ftsdf(k+1) = -euler(4,n2)*distd(3,n) + euler(7,n2)*distd(2,n)
       ftsdf(k+2) = -euler(7,n2)*distd(1,n) + euler(1,n2)*distd(3,n)
       ftsdf(k+3) = -euler(1,n2)*distd(2,n) + euler(4,n2)*distd(1,n)

       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(2,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1 + ndime   !compute next pointer
       nesdf(k)   = 2+10*n2         !store master node and dof until later
       nesdf(k+1) = 4+10*n2         !
       nesdf(k+2) = 5+10*n2         !
       nesdf(k+3) = 6+10*n2         !
       ftsdf(k)   = 1d0             !factor for translational DOF & rotational DOFs
       ftsdf(k+1) = -euler(5,n2)*distd(3,n) + euler(8,n2)*distd(2,n)
       ftsdf(k+2) = -euler(8,n2)*distd(1,n) + euler(2,n2)*distd(3,n)
       ftsdf(k+3) = -euler(2,n2)*distd(2,n) + euler(5,n2)*distd(1,n)

       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(3,n1)= -j              !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1 + ndime   !compute next pointer
       nesdf(k)   = 3+10*n2         !store master node and dof until later
       nesdf(k+1) = 4+10*n2         !
       nesdf(k+2) = 5+10*n2         !
       nesdf(k+3) = 6+10*n2         !
       ftsdf(k)   = 1d0             !factor for translational DOF & rotational DOFs
       ftsdf(k+1) = -euler(6,n2)*distd(3,n) + euler(9,n2)*distd(2,n)
       ftsdf(k+2) = -euler(9,n2)*distd(1,n) + euler(3,n2)*distd(3,n)
       ftsdf(k+3) = -euler(3,n2)*distd(2,n) + euler(6,n2)*distd(1,n)
       ! IF slave node have local system
       IF( naeul(n1) )THEN   ! rotational degrees of freedom
         k = npsdf(j)               !first pointer
         CALL proma4(ftsdf(k),euler(1,n2),euler(1,n1),3,3,3)
         DO i=1,ndime                 !for each slave rotational DOF
           k = npsdf(j)               !first pointer
           ifpre(ndime+i,n1) = -j     !slave DOF order
           j = j+1                    !increase counter of slave DOFs
           npsdf(j)   = k + ndime     !compute next pointer
           DO l=1,ndime      !for each master rotational DOF
             nesdf(k) = ndime+l+10*n2 !store master node and dof until later
             k = k+1                  !counter in NESDF & FTSDF
           END DO
         END DO
       ELSE
         ifpre(4:ndofn,n1) = 1              !no DOF
       END IF
     ELSE ! ndime = 2 (2-D problems)
       t1 = COS(euler(1,n2))                !local system of master node
       t2 = SIN(euler(1,n2))
       distd(1,n)= t1*d(1)+t2*d(2)          !material components
       distd(2,n)=-t2*d(1)+t1*d(2)          !of distance vector
       k = npsdf(j)                         !pointer to first slave DOFs
       ifpre(1,n1) = -j                     !slave DOF order
       j = j+1                              !increase counter of slave DOFs
       npsdf(j)   = k + 2                   !compute next pointer
       nesdf(k)   = 1+10*n2                 !store master node and dof until later
       nesdf(k+1) = 3+10*n2                 !store master node and dof until later
       ftsdf(k)   = 1d0                     !factor for translational DOF
       ftsdf(k+1) = -t2*distd(1,n) - t1*distd(2,n) !factor for rotational DOF

       k = npsdf(j)                         !pointer to first slave DOFs
       ifpre(2,n1) = -j                     !slave DOF order
       j = j+1                              !increase counter of slave DOFs
       npsdf(j)   = k + 2                   !compute next pointer
       nesdf(k)   = 2+10*n2                 !store master node and dof until later
       nesdf(k+1) = 3+10*n2                 !store master node and dof until later
       ftsdf(k)   = 1d0                     !factor for translational DOF
       ftsdf(k+1) = t1*distd(1,n) - t2*distd(2,n) !factor for rotational DOF
       ! IF slave node have local system
       IF( naeul(n1) )THEN                  ! rotational degree of freedom
         distd(3,n)= euler(1,n1)-euler(1,n2)  !relative angle
         k = npsdf(j)                         !pointer to first slave DOFs
         ifpre(3,n1) = -j                     !slave DOF order
         j = j+1                              !increase counter of slave DOFs
         npsdf(j)   = k + 1                   !compute next pointer
         nesdf(k)   = 3+10*n2                 !store master node and dof until later
         ftsdf(k)   = 1d0                     !factor for rotational DOF
       ELSE
         ifpre(3:ndofn,n1) = 1              !no DOF
       END IF
     END IF
   ELSE  !model does not include rotational DOFs
     distd(1:ndime,n) = d(1:ndime) !material components
     ! translational degrees of freedom
     DO i=1,ndime
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(i,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1           !compute next pointer
       nesdf(k)   = i+10*n2         !store master node and dof until later
       ftsdf(k)   = 1d0             !factor for translational DOF
     END DO
   END IF

 END DO
 IF( nndp == ndepd )RETURN
 ! rigid body nodes
 n = nndp              !last position
 DO n1=1,npoin         !loop for each dependant node
   IF( ifpre(1,n1) >= -nn )CYCLE
   n  = n + 1             !increase counter
   n2 = -ifpre(1,n1)-nn   !master node label
   nndpd(1,n) = n1              !copy slave into array
   nndpd(2,n) = n2              !copy master into array
   DO i = 1,nrotd               !for each (slave) DOF
     IF(ifpre(i,n1) >= -nn ) THEN  !check that the DOF is not already slave
       WRITE(lures,"(' DOF already declared Slave:  ',&
                  &  'Node',i6,' DOF',i2)",ERR=9999) label(n1) ,i
       CALL runend('READPD: Inconsistent Input Data    ')
     END IF
   END DO
   ! factors of commented lines depends on coordinates of the master node (not computed yet)
   IF( eul )THEN
     IF(ndime == 3) THEN            !for 3-D problems
       !distd(1:3,n) = MATMUL(d,RESHAPE(euler(1:9,n2),(/3,3/))) !material components
       ! translational degrees of freedom
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(1,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1 + ndime   !compute next pointer
       nesdf(k)   = 1+10*n2         !store master node and dof until later
       nesdf(k+1) = 4+10*n2         !
       nesdf(k+2) = 5+10*n2         !
       nesdf(k+3) = 6+10*n2         !
       ftsdf(k)   = 1d0             !factor for translational DOF & rotational DOFs
       !ftsdf(k+1) = -euler(4,n2)*distd(3,n) + euler(7,n2)*distd(2,n)
       !ftsdf(k+2) = -euler(7,n2)*distd(1,n) + euler(1,n2)*distd(3,n)
       !ftsdf(k+3) = -euler(1,n2)*distd(2,n) + euler(4,n2)*distd(1,n)

       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(2,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1 + ndime   !compute next pointer
       nesdf(k)   = 2+10*n2         !store master node and dof until later
       nesdf(k+1) = 4+10*n2         !
       nesdf(k+2) = 5+10*n2         !
       nesdf(k+3) = 6+10*n2         !
       ftsdf(k)   = 1d0             !factor for translational DOF & rotational DOFs
       !ftsdf(k+1) = -euler(5,n2)*distd(3,n) + euler(8,n2)*distd(2,n)
       !ftsdf(k+2) = -euler(8,n2)*distd(1,n) + euler(2,n2)*distd(3,n)
       !ftsdf(k+3) = -euler(2,n2)*distd(2,n) + euler(5,n2)*distd(1,n)

       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(3,n1)= -j              !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1 + ndime   !compute next pointer
       nesdf(k)   = 3+10*n2         !store master node and dof until later
       nesdf(k+1) = 4+10*n2         !
       nesdf(k+2) = 5+10*n2         !
       nesdf(k+3) = 6+10*n2         !
       ftsdf(k)   = 1d0             !factor for translational DOF & rotational DOFs
       !ftsdf(k+1) = -euler(6,n2)*distd(3,n) + euler(9,n2)*distd(2,n)
       !ftsdf(k+2) = -euler(9,n2)*distd(1,n) + euler(3,n2)*distd(3,n)
       !ftsdf(k+3) = -euler(3,n2)*distd(2,n) + euler(6,n2)*distd(1,n)
       ! IF slave node have local system
       IF( naeul(n1) )THEN   ! rotational degrees of freedom
         DO i=1,ndime                 !for each rotational DOF
           m=3*i-2                    !first position in rotation matrix
           k = npsdf(j)               !first pointer
           ifpre(ndime+i,n1) = -j     !slave DOF order
           j = j+1                    !increase counter of slave DOFs
           npsdf(j)   = k + ndime     !compute next pointer
           DO l=1,ndime
             ll = 3*l-2
             nesdf(k) = ndime+l+10*n2 !store master node and dof until later
             !l_ij = lM_i . lS_j
             ftsdf(k) = DOT_PRODUCT(euler(m:m+2,n1),euler(ll:ll+2,n2)) !factor
             k = k+1                  !counter in NESDF & FTSDF
           END DO
         END DO
       ELSE
         ifpre(4:ndofn,n1) = 1              !no DOF
       END IF
     ELSE ! ndime = 2 (2-D problems)
       !t1 = COS(euler(1,n2))                !local system of master node
       !t2 = SIN(euler(1,n2))
       !distd(1,n)= t1*d(1)+t2*d(2)          !material components
       !distd(2,n)=-t2*d(1)+t1*d(2)          !of distance vector
       k = npsdf(j)                         !pointer to first slave DOFs
       ifpre(1,n1) = -j                     !slave DOF order
       j = j+1                              !increase counter of slave DOFs
       npsdf(j)   = k + 2                   !compute next pointer
       nesdf(k)   = 1+10*n2                 !store master node and dof until later
       nesdf(k+1) = 3+10*n2                 !store master node and dof until later
       ftsdf(k)   = 1d0                     !factor for translational DOF
       !ftsdf(k+1) = -t2*distd(1,n) - t1*distd(2,n) !factor for rotational DOF

       k = npsdf(j)                         !pointer to first slave DOFs
       ifpre(2,n1) = -j                     !slave DOF order
       j = j+1                              !increase counter of slave DOFs
       npsdf(j)   = k + 2                   !compute next pointer
       nesdf(k)   = 2+10*n2                 !store master node and dof until later
       nesdf(k+1) = 3+10*n2                 !store master node and dof until later
       ftsdf(k)   = 1d0                     !factor for translational DOF
       !ftsdf(k+1) = t1*distd(1,n) - t2*distd(2,n) !factor for rotational DOF
       ! IF slave node have local system
       IF( naeul(n1) )THEN                  ! rotational degree of freedom
         !distd(3,n)= euler(1,n1)-euler(1,n2)  !relative angle
         k = npsdf(j)                         !pointer to first slave DOFs
         ifpre(3,n1) = -j                     !slave DOF order
         j = j+1                              !increase counter of slave DOFs
         npsdf(j)   = k + 1                   !compute next pointer
         nesdf(k)   = 3+10*n2                 !store master node and dof until later
         ftsdf(k)   = 1d0                     !factor for rotational DOF
       ELSE
         ifpre(3,n1) = 1              !no DOF
       END IF
     END IF
   ELSE  !model does not include rotational DOFs
     ! distd(1:ndime,n) = d(1:ndime) !material components
     ! translational degrees of freedom
     DO i=1,ndime
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(i,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1           !compute next pointer
       nesdf(k)   = i+10*n2         !store master node and dof until later
       ftsdf(k)   = 1d0             !factor for translational DOF
     END DO
     ifpre(ndime+1:nrotd,n1) = 1    !No DOF
   END IF
   IF( n == ndepd)EXIT
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE readpd
