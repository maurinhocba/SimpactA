 SUBROUTINE arisfi(naris,x,nndpd,euler)
 !***********************************************************************
 !
 !     calculates dependant displacements
 !
 !***********************************************************************
 USE npo_db,  ONLY : naeul,ifpre
 USE ctrl_db, ONLY : ndime,neulr
 USE kinc_db, ONLY : npsdf,nesdf,ftsdf,nn
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: naris,nndpd(3,naris)
 REAL (kind=8),INTENT(IN OUT) :: x(:,:)
 REAL (kind=8),POINTER :: euler(:,:)

 REAL(kind=8), PARAMETER :: minp = 0.5d0            !minimun proyection allowed (60 degrees)
 INTEGER(kind=4), PARAMETER :: ll(2,3) = (/ 2,3,  1,3,  1,2 /)  !slave DOFs
 INTEGER (kind=4) :: i,j,k,l,ii,kk,n,ns,n1,n2,m,mm,jj,lp(2),idof
 REAL (kind=8) :: d(ndime),lg,phi(3),alpha,t1
 LOGICAL :: eul,seul,slide

 INTERFACE
   INCLUDE 'actrot.h'
 END INTERFACE

 eul = neulr > 0  !local system exist
 seul = .FALSE.   !initializes
 DO m=1,naris                 !for each node on a side
   ns = nndpd(1,m)            !slave node
   n1 = nndpd(2,m)            !first master node
   n2 = nndpd(3,m)            !second master node
   slide = ns < 0
   IF(slide) ns = ABS(ns)
   d = x(1:ndime,n2) - x(1:ndime,n1)  !present side
   CALL vecuni(ndime,d,lg)     !unit side vector
   IF( slide )THEN
     alpha = DOT_PRODUCT(x(:,ns)-x(:,n1),d)/lg
   ELSE
     k = npsdf(-ifpre(1,ns))
     alpha = (1d0-ftsdf(k))
   END IF
   x(1:ndime,ns) = x(1:ndime,n1) + d * lg * alpha  !present slave node position
   IF( slide )THEN       !recompute tangent factors

     DO i=1,ndime      !search present free DOF as it is not in data base
       IF( ifpre(i,ns) > 0  .OR. ifpre(i,ns) < -nn )THEN   ! free DOFs
         mm = i         ! free DOF
         EXIT          ! once found exit loop
       END IF
     END DO
     idof= ifpre(mm,ns)              !keep internal DOF
     IF( idof > 0 )THEN             !only for free DOF
       IF( ABS(d(mm))  < minp )THEN !change DOF if necessary
         DO l=1,ndime-1
           lp(l) = ifpre(ll(l,mm),ns)     !pointers to slave DOFs
           IF( ABS(d(ll(l,mm))) > ABS(d(mm)) ) j=ll(l,mm)
         END DO
         DO jj=1,ndime-1   !for each slave DOF
           ii = ll(jj,j)                 !new slave DOF
           ifpre(ii,n1) = lp(jj)         !pass pointer to npsdf
           k = npsdf(-lp(jj))            !position in npsdf of old dependency
           nesdf(k)   = ifpre(ii,n1)     !recompute master DOFs
           nesdf(k+1) = ifpre(ii,n2)     !
           !nesdf(k+2) = idof             !does not change
           nesdf(k+3) = ifpre(j,n1)      !
           nesdf(k+4) = ifpre(j,n2)      !
         END DO
         mm = j                        !permute free DOF
         ifpre(mm,ns) = idof           !pass internal DOF
       END IF
     END IF

     lp = ll(:,mm)   ! slave DOFs
     DO i=1,ndime-1                !for each translational slave DOF
       t1 = d(lp(i)) / d(mm)          !ti/lm
       j = - ifpre(lp(i),ns)         !slave DOF order (position)
       k = npsdf(j)                  !first position for this DOF in NESDF & FTSDF
       ftsdf(k)   = 1d0-alpha        !factor for v i on master(1st node)
       ftsdf(k+1) = alpha            !factor for v i on master(2nd node)
       ftsdf(k+2) = t1               !factor for v m on slave
       ftsdf(k+3) = -t1*(1d0-alpha)  !factor for v i on master(1st node)
       ftsdf(k+4) = -t1*alpha        !factor for v i on master(2nd node)
     END DO

   ELSE
     IF( eul ) seul = naeul(ns) !see if local system exists at slave node
     IF( .NOT.seul )CYCLE       !cycle if no local system

     ! update rotational DOFs dependency factors
     j = -ifpre(ndime+1,ns)     !rotational dof's position
     IF(ndime == 2) THEN
       !  computes the new angle
       IF(d(1) <= 0)euler(1,ns) =  ACOS(d(2))
       IF(d(1) > 0 )euler(1,ns) = -ACOS(d(2))
       k = npsdf(j)
       ftsdf(k)   =  d(2)*l
       ftsdf(k+1) =  d(1)*l
       ftsdf(k+2) = -d(2)*l
       ftsdf(k+3) = -d(1)*l
       j = j+1
     ELSE
       !  computes the new cartesian system
       CALL vecpro(euler(7:9,ns),d,phi)
       CALL proma1(d,phi,euler(1,ns),1,3,3)
       CALL actrot(euler(1:9,ns),d(1:3))
       !   generates  factors
       jj = npsdf(j)
       kk = 3
       DO k=4,5
         DO n=1,2
           DO i=1,ndime
             ftsdf(jj) = l*euler(kk+i,ns)
             jj = jj+1
           END DO
           l = -l
         END DO
         kk = 0
         l = -l
       END DO
       j = j+2
     END IF
   END IF
 END DO
 RETURN

 END SUBROUTINE arisfi
