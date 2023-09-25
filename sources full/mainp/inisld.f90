 SUBROUTINE inisld(npsdf,nesdf,ftsdf,j,np)
 !**********************************************************************
 !
 !     APPLIES and generates DATA for sliding nodes
 !     for 3-D problems
 !
 !**********************************************************************
 USE lispa0, ONLY : lures
 USE nsld_db
 USE ctrl_db, ONLY : ndime,ndofn,npoin,neulr
 USE npo_db,  ONLY : naeul,ifpre,coora,euler,label,eule0
 USE kinc_db, ONLY : sldnd,sldvl
 USE outp_db, ONLY: iwrit

 IMPLICIT NONE
 ! dummy arguments
 REAL(kind=8), PARAMETER :: toler = 1.d0
 INTEGER (kind=4),INTENT(IN) :: np    !size (to check it is not exceeded)
 INTEGER (kind=4),INTENT(IN OUT) :: j,npsdf(np),nesdf(4*np)
 REAL (kind=8), INTENT(IN OUT) :: ftsdf(4*np)
 ! local variables
 INTEGER (kind=4) :: i,k,l(2),m,n,n1,n2
 REAL (kind=8) :: d(ndime),t1,t2,ang(3),t(3)
 LOGICAL :: sld,rot

 INTEGER(kind=4) chnode ! function
 INTERFACE
   INCLUDE 'inrotm.h'
  END INTERFACE

 IF( (ndime == 3 .AND. ndofn < 6 ).AND. iwrit == 1) &
    WRITE(lures,"(' its not posible to have true sliding nodes in',&
          &    ' 3-D with only',i2,' DOF')",ERR=9999) ndofn

 IF (iwrit == 1 .AND. nnsld > 0 ) WRITE(lures,"(//,3x,'Slave/Master nodes Data',//&
                 & '   Slave Node     Master Node')",ERR=9999)

 DO n=1,nnsld          !loop for each dependant node

   n1 =nsldat(1,n)    !slave node label
   n2 =nsldat(2,n)    !master node label
   sld = .NOT.nsldfg(2,n)   !sliding
   rot = .NOT.nsldfg(3,n)   !rotating
   IF (iwrit == 1) WRITE(lures, '(2i10)',ERR=9999) n1, n2   !ECHO
   n1 = chnode(n1)              !slave internal node
   n2 = chnode(n2)              !master internal node
   sldnd(1,n) = n1              !copy slave into array
   sldnd(2,n) = n2              !copy master into array
   sldnd(3,n) = 1                               !all flags are FALSE
   IF( nsldfg(1,n) )sldnd(3,n) = -1             !AUTOM flag is TRUE
   IF( nsldfg(2,n) )sldnd(3,n) = sldnd(3,n) * 2 !NOSLD flag is TRUE
   IF( nsldfg(3,n) )sldnd(3,n) = sldnd(3,n) * 3 !NOROT flag is TRUE
   DO i = 1,ndofn               !for each (slave) DOF
     IF(ifpre(i,n1) < 0 ) THEN  !check that the DOF is not already slave
       WRITE(lures,"(' DOF already declared Slave:  ',&
                  &  'Node',i6,' DOF',i2)",ERR=9999) label(n1) ,i
       CALL runend('READPD: Inconsistent Input Data    ')
     END IF
   END DO
   d = coora(1:ndime,n1) - coora(1:ndime,n2)    !distance vector between nodes (master->slave)
   IF( .NOT. naeul(n2) )THEN
     naeul(n2) = .TRUE.             !local system required at master node
     IF( ALL(euler(4:9,n2) == 0d0))THEN
       ang = euler(1:3,n2)            !Recover the Euler angles (radians)
       CALL inrotm(ang,euler(:,n2))   !compute rotation matrix
     END IF
   END IF
   IF( .NOT. naeul(n1) )THEN          !If no local system at slave node
     IF( .NOT. nsldfg(3,n) )THEN      !If relative rotations are considered
       naeul(n1) = .TRUE.             !local system required at slave node
       IF( ALL(euler(4:9,n1) == 0d0))THEN
         ang = euler(1:3,n1)              !Recover the Euler angles (radians)
         CALL inrotm(ang,euler(:,n1))     !compute rotation matrix
       END IF
     END IF
   END IF

   ! checks

   IF( naeul(n1) ) THEN
     !verify T1M == T1S
     t1 = DOT_PRODUCT( euler(1:3,n1),euler(1:3,n2) )  !cos T1M T1S
     IF( t1 < toler ) THEN                     !not equal 1
       WRITE(*,*)' directions 1 in nodes ',n1,n2,' NOT parallel, MODIFIED'
       !modify
       CALL vecpro(euler(1,n1),euler(1,n2),t(1))    !normal to both t1 vectors
       CALL vecuni(3,t(1),t1)                       !angle between then
       t1 = ASIN(t1)
       ang(1) = 0d0                                 !rotation vector between systems
       ang(2) = t1*DOT_PRODUCT(t,euler(4:6,n1))
       ang(3) = t1*DOT_PRODUCT(t,euler(7:9,n1))
       CALL actrot(euler(1,n1),ang(1))              !update local system at slave node
       CALL angeul(euler(1:9,n1),eule0(1:3,n1),.TRUE.) !regenerates Euler angles (in Rads)
     END IF
     t1 = DOT_PRODUCT(euler(4:6,n1),euler(4:6,n2))     !cos
     t2 = DOT_PRODUCT(euler(4:6,n1),euler(7:9,n2))     !sin
     sldvl(2,n) = ATAN2(t2,t1)                         !keep angle
   ELSE
     sldvl(2,n) = 0                                    !no angle
   END IF
   !check alignment
   t1 = SQRT(DOT_PRODUCT(d,d))    !distance length
   IF( t1 > 0d0 )THEN             !if nodes not coincident
     t2 = DOT_PRODUCT(d,euler(1:3,n2))     !proyection
     !check parallelism
     coora(:,n1) = coora(:,n2) + t2*euler(1:3,n2)     !modify original coordinates of slave
     IF( ABS(t2/t1) < toler ) THEN
       WRITE(*,*)' nodes slave and master ',n1,n2,' NOT aligned, MODIFIED'
       WRITE(lures,*)' nodes slave and master ',n1,n2,' NOT aligned, MODIFIED'
     END IF       
     WRITE(lures,"(i6,3e24.16)")n1,coora(:,n1)
   ELSE
     t2 = 0d0
   END IF
   sldvl(1,n) = t2              !distance

   IF( sld )THEN   ! sliding
     m = 1 ; l =(/2,3/)        ! default free DOF and slave DOFs
     IF( nsldfg(1,n) )THEN      !IF AUTOMatic
       IF( ABS(euler(2,n2)) > ABS(euler(1,n2)) )THEN !compare proyections along x1 and x2
          m=2 ; l(1) = 1                 !change if x2 > x1
       END IF
       IF( ABS(euler(3,n2)) > ABS(euler(m,n2)) )THEN !compare proyection along x3
          m=3 ; l=(/1,2/)                !change if x3 is the largent
       END IF
     END IF
     DO i=1,2 !for each slave translational DOF
       t1 = euler(l(i),n2)/euler(m,n2) !li1/lm1
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(l(i),n1) = -j          !store slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 5           !compute next pointer
       nesdf(k)   = m+10*n1         !store master node and dof until later
       nesdf(k+1) = m+10*n2         !
       nesdf(k+2) = l(i)+10*n2      !
       nesdf(k+3) = 5+10*n2         !
       nesdf(k+4) = 6+10*n2         !
       ftsdf(k)   =  t1                       !factor for v m on slave
       ftsdf(k+1) = -t1                       !factor for v m on master
       ftsdf(k+2) = 1d0                       !factor for v l on master
       ftsdf(k+3) = (-euler(6+l(i),n2) + euler(6+m,n2)*t1 )*t2  !factor for w2 on master
       ftsdf(k+4) = ( euler(3+l(i),n2) - euler(3+m,n2)*t1 )*t2  !factor for w3 on master
     END DO

   ELSE         !no sliding (but rotating)
     DO i=1,3                   ! for each translational DOF
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(i,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 3           !compute next pointer
       nesdf(k)   = i+10*n2         !store master node and dof until later
       nesdf(k+1) = 5+10*n2         !
       nesdf(k+2) = 6+10*n2         !
       ftsdf(k)   = 1d0                   !factor for v_i on master
       ftsdf(k+1) = -euler(6+i,n2)*t2     !factor for w2 on master
       ftsdf(k+2) =  euler(3+i,n2)*t2     !factor for w3 on master
     END DO
   END IF

   IF( rot )THEN   !rotating allowed
     t1 = DOT_PRODUCT(euler(4:6,n1),euler(4:6,n2)) !cos phi
     t2 = DOT_PRODUCT(euler(4:6,n1),euler(7:9,n2)) !sin phi
     DO i=5,6      !for each slave rotational DOF
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(i,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 2           !compute next pointer
       nesdf(k)   = 5+10*n2         !store master node and dof until later
       nesdf(k+1) = 6+10*n2         !
       ftsdf(k)   =  t1                       !factor for w2 on slave  (COS)
       ftsdf(k+1) =  t2                       !factor for w3 on master (SIN)
       t1 = -t2                     !rotate 90 degrees for next DOF
       t2 = ftsdf(k)
     END DO

   ELSE IF( naeul(n1) )THEN         !no relative rotation allowed
     DO i=4,6                  !for each rotational DOF
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(i,n1) = -j             !slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 1           !compute next pointer
       nesdf(k)   = i+10*n2         !store master node and dof until later
       ftsdf(k)   =  1d0            !factor for w1 on slave
     END DO
   END IF
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inisld
