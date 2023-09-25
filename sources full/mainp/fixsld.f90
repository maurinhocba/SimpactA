 SUBROUTINE fixsld(x,euler)
 !***********************************************************************
 !
 !     updates configuration and factors of sliding node
 !
 !***********************************************************************
 USE kinc_db, ONLY : sldnd,sldvl,npsdf,nesdf,ftsdf,distd,nn
 USE npo_db,  ONLY : naeul,ifpre
 USE nsld_db, ONLY : nnsld
 IMPLICIT NONE
 ! dummy arguments
 REAL (kind=8),INTENT(IN OUT) :: x(:,:)  !present coordinates
 REAL (kind=8), POINTER :: euler(:,:)    !present local system
 ! local variables
 REAL(kind=8), PARAMETER :: toler = 0.9999999d0,  & !to compare with 1.0
                            minp = 0.5d0            !minimun proyection allowed (60 degrees)
 INTEGER(kind=4), PARAMETER :: ll(2,3) = (/ 2,3,  1,3,  1,2 /)  !slave DOFs
 INTEGER (kind=4) :: i,j,k,n,n1,n2,m,l(2),ii,jj,idof
 REAL    (kind=8) :: d(3),t1,t2,ang(3),dd
 LOGICAL :: aut,sld,rot

 DO n=1,nnsld       !for each sliding node
   n1 = sldnd(1,n)     !slave node
   n2 = sldnd(2,n)     !master node
   aut = sldnd(3,n) < 0        !automatic selection of free translational DOF
   sld = ABS(sldnd(3,n)) /= 2  !sliding
   rot = ABS(sldnd(3,n)) /= 3  !rotating
   !recompute slave node local system
   IF( naeul(n1) .AND. rot ) THEN ! if relative rotations are possible
     t1 = DOT_PRODUCT( euler(1:3,n1),euler(1:3,n2) )  !angle cosine between T1 axes
     IF( t1 < toler ) THEN                            !if less than 1, i.e. not parallel
       !modify
       CALL vecpro(euler(1,n1),euler(1,n2),d(1))    !normal to both t1 vectors
       CALL vecuni(3,d(1),t1)                       !angle between then
       t1 = ASIN(t1)
       ang(1) = 0d0                                 !rotation vector components
       ang(2) = t1*DOT_PRODUCT(d,euler(4:6,n1))
       ang(3) = t1*DOT_PRODUCT(d,euler(7:9,n1))
       CALL actrot(euler(1,n1),ang(1))              !modify slave local system to match Master t1
     END IF
     t1 = DOT_PRODUCT(euler(4:6,n1),euler(4:6,n2))     !cos
     t2 = DOT_PRODUCT(euler(4:6,n1),euler(7:9,n2))     !sin
     sldvl(2,n) = ATAN2(t2,t1)                         !keep angle (for reference, not used)
   END IF
   !recompute slave node position
   IF( sld ) THEN ! sliding
     d = x(1:3,n1) - x(1:3,n2)             !distance vector between nodes (master->slave)
     dd = DOT_PRODUCT(d,euler(1:3,n2))     !proyection
     sldvl(1,n) = dd                       !distance
   ELSE           ! no sliding
     dd = sldvl(1,n)    !recove constant distance
   END IF
   x(:,n1) = x(:,n2) + dd*euler(1:3,n2)     !modify original coordinates of slave

   !recompute tangent factors
   IF( sld  )THEN        ! sliding
     IF( aut )THEN   !automatic
       DO i=1,3          !search present free DOF as it is not in data base
         IF( ifpre(i,n1) > 0  .OR. ifpre(i,n1) < -nn )THEN   ! free DOFs
           m = i         ! free DOF
           EXIT          ! once found exit loop
         END IF
       END DO
       idof= ifpre(m,n1)              !keep internal DOF
       IF( idof > 0 )THEN             !only for free DOF
         IF( ABS(euler(m,n2)) < minp )THEN !change DOF if necessary
           l = ifpre(ll(:,m),n1)     !pointers to slave DOFs
           IF( ABS(euler(ll(1,m),n2)) > ABS(euler(      m,n2)) ) j=ll(1,m)
           IF( ABS(euler(ll(2,m),n2)) > ABS(euler(ll(1,m),n2)) ) j=ll(2,m)
           DO jj=1,2   !for each slave DOF
             ii = ll(jj,j)                 !new slave DOF
             ifpre(ii,n1) = l(jj)          !pass pointer to npsdf
             k = npsdf(-l(jj))             !position in npsdf of old dependency
             nesdf(k)   = idof             !recompute master DOFs
             nesdf(k+1) = ifpre(j,n2)      !
             nesdf(k+2) = ifpre(ii,n2)     !
           END DO
           m = j                        !permute free DOF
           ifpre(m,n1) = idof           !pass internal DOF
         END IF
       END IF
     ELSE
       m = 1     ! default free DOFs
     END IF
     l = ll(:,m)   ! slave DOFs
     DO i=1,2                 !for each translational slave DOF
       j = - ifpre(l(i),n1)            !slave DOF order (position)
       k = npsdf(j)                    !first position for this DOF in NESDF & FTSDF
       t1 = euler(l(i),n2)/euler(m,n2)       !ll1/lm1
       ftsdf(k)   =  t1                       !factor for v m on slave
       ftsdf(k+1) = -t1                       !factor for v m on master
       ftsdf(k+3) = (-euler(6+l(i),n2) + euler(6+m,n2)*t1 )*dd  !factor for w2 on master
       ftsdf(k+4) = ( euler(3+l(i),n2) - euler(3+m,n2)*t1 )*dd  !factor for w3 on master
     END DO
   ELSE                 ! no sliding
     DO i=1,3  !for each translational DOF
       j = - ifpre(i,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k+1) = -euler(6+i,n2)*dd     !factor for w2 on master
       ftsdf(k+2) =  euler(3+i,n2)*dd     !factor for w3 on master
     END DO
   END IF

   IF( rot ) THEN ! rotating
     t1 = DOT_PRODUCT(euler(4:6,n1),euler(4:6,n2)) !cos phi
     t2 = DOT_PRODUCT(euler(4:6,n1),euler(7:9,n2)) !sin phi
     DO i=5,6
       j = - ifpre(i,n1)            !slave DOF order
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ftsdf(k)   =  t1                       !factor for w2 on slave
       ftsdf(k+1) =  t2                       !factor for w3 on master
       t1 = -t2                     !rotate 90 degrees
       t2 = ftsdf(k)
     END DO
   END IF

 END DO
 RETURN

 END SUBROUTINE fixsld
