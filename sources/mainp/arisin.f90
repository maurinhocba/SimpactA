 SUBROUTINE arisin(naris,npsdf,nesdf,ftsdf,nndpd,j,np)
 !**********************************************************************
 !
 !     generates DATA for nodes fixed on a side
 !
 !**********************************************************************
 USE nar_db
 USE lispa0
 USE ctrl_db, ONLY : ndime,nrotd,npoin,neulr
 USE npo_db, ONLY : ifpre,coora,euler,label,naeul
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: naris,np
 INTEGER (kind=4),INTENT(IN OUT) :: npsdf(np),j
 INTEGER (kind=4),INTENT(OUT) :: nesdf(2*np),nndpd(3,naris)
 REAL    (kind=8),INTENT(OUT) :: ftsdf(2*np)

 INTEGER (kind=4) :: i,jj,j0,k,kk,m,n,nn(2),ns,mm,ll(2)
 REAL    (kind=8) :: alpha,l,d(ndime),phi(3),t1

 INTEGER(kind=4) chnode
 LOGICAL :: eul,slide

 INTERFACE
   INCLUDE 'actrot.h'
 END INTERFACE

 WRITE(lures,"(' Slave Nodes Data',/,'  Slave Node   Side Nodes')",ERR=9999)

 eul = neulr > 0
 DO m=1,naris                 !for each slave NODE
   ns   =nardat(1,m)          !slave node label
   nn(1)=nardat(2,m)          !first master node label
   nn(2)=nardat(3,m)          !second master node label
   slide = ns < 0
   WRITE(lures,"(i12,2i6)",ERR=9999) ABS(ns),nn(1:2)
   nn(1) = chnode(nn(1))  !first master node internal number
   nn(2) = chnode(nn(2))  !second master node internal number
   ns    = chnode(ABS(ns))  !slave node internal number
   DO i = 1,nrotd             !Translational DOFs
     IF(ifpre(i,ns) < 0) THEN
       WRITE(lures,"(' slave DOF already declared Node',i6,' dof',i2)", &
             ERR=9999) label(ns) ,i
       CALL runen3('ARISIN: Inconsistent Input Data    ')
     END IF
   END DO
   nndpd(1,m) = ns                    !store slave node in array
   IF(slide) nndpd(1,m) = -ns         !flag
   nndpd(2:3,m) = nn(1:2)             !store master nodes in array
   !  proyects the node on the side
   d = coora(1:ndime,nn(2)) -  coora(1:ndime,nn(1))   !arista n1 --> n2
   l = DOT_PRODUCT(d,d)     ! squared lenght of the arista
   ! relative position from n1
   alpha = DOT_PRODUCT(d,coora(1:ndime,ns)-coora(1:ndime,nn(1)))/l
   l    = 1d0/SQRT(l)                              !side lenght  ^(-1)
   !  correct coordinates of slave node to be on the side
   coora(1:ndime,ns) = coora(1:ndime,nn(1)) + alpha*d
   d = d*l     !unit side vector
   IF( slide ) THEN             !if node can slide along side
     IF( ndime == 3 )THEN
       mm = 3 ; ll =(/1,2/)        ! default free DOF and slave DOFs
       IF( ABS(d(2)) > ABS(d(3)) )THEN !compare proyection along x2
          mm=2 ; ll(2)=3                !change if x2 is the largent
       END IF
       IF( ABS(d(1)) >  ABS( d(mm) ) )THEN !compare proyections along x1
          mm=1 ; ll = (/2,3/)              !change if x1 > x3
       END IF

       DO i=1,2 !for each slave translational DOF
         t1 = d(ll(i)) / d(mm)        !ti/lm
         k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
         ifpre(ll(i),ns) = -j         !store slave DOF order
         j = j+1                      !increase counter of slave DOFs
         npsdf(j)   = k + 5           !compute next pointer
         nesdf(k)   = ll(i)+10*nn(1)  !store master node and dof until later
         nesdf(k+1) = ll(i)+10*nn(2)  !
         nesdf(k+2) = mm   +10*ns     !
         nesdf(k+3) = mm   +10*nn(1)  !
         nesdf(k+4) = mm   +10*nn(2)  !
         ftsdf(k)   = 1d0-alpha        !factor for v i on master(1st node)
         ftsdf(k+1) = alpha            !factor for v i on master(2nd node)
         ftsdf(k+2) = t1              !factor for v m on slave
         ftsdf(k+3) = -t1*(1d0-alpha)  !factor for v i on master(1st node)
         ftsdf(k+4) = -t1*alpha        !factor for v i on master(2nd node)
       END DO
     ELSE
       mm= 2 ; ll(1) = 1        ! default free DOF and slave DOFs
       IF( ABS(d(1)) > ABS(d(2)) )THEN !compare proyection along x1
          mm=1 ; ll(1)=2               !change if x1 is the largent
       END IF

       t1 = d(ll(1)) / d(mm)        !ti/lm
       k = npsdf(j)                 !first position for this DOF in NESDF & FTSDF
       ifpre(ll(i),nn(1)) = -j         !store slave DOF order
       j = j+1                      !increase counter of slave DOFs
       npsdf(j)   = k + 5           !compute next pointer
       nesdf(k)   = ll(1)+10*nn(1)  !store master node and dof until later
       nesdf(k+1) = ll(1)+10*nn(2)  !
       nesdf(k+2) = mm   +10*ns     !
       nesdf(k+3) = mm   +10*nn(1)  !
       nesdf(k+4) = mm   +10*nn(2)  !
       ftsdf(k)   = 1d0-alpha        !factor for v i on master(1st node)
       ftsdf(k+1) = alpha            !factor for v i on master(2nd node)
       ftsdf(k+2) = t1               !factor for v m on slave
       ftsdf(k+3) = -t1*(1d0-alpha)  !factor for v i on master(1st node)
       ftsdf(k+4) = -t1*alpha        !factor for v i on master(2nd node)

     END IF
   ELSE
     !  generates connections and factors for trans. DOF
     DO i = 1,ndime             !Translational DOFs
       jj = npsdf(j)                 !initial position in long vector
       ifpre(i,ns) = -j
       DO n=1,2                      !for each master node
         alpha = 1d0-alpha          !permutes factor
         nesdf(jj) = i+10*nn(n)      !store master node and dof until later
         ftsdf(jj) = alpha            !factor
         jj = jj+1                   !next position in long vector
       END DO
       j = j + 1              !next slave DOF
       IF( j > np )CALL runen3('ARISIN: insufficient auxiliar space (MEM)')
       npsdf(j) = jj          !pointer of next slave DOF to long vector
     END DO
     !  generates connections and factors for Rotational DOFs
     IF( eul )THEN          !if Local systems exist
       IF( naeul(ns) )THEN  !if the node have local system
         DO k = 1,ndime-1            !Rotational DOFs
           i = k+ndime
           jj = npsdf(j)               !initial position in long vector
           ifpre(i,ns) = -j
           DO n=1,2                      !for each master node
             DO i=1,ndime                !for each direction
               nesdf(jj) = i+10*nn(n)    !store master node and dof until later
               jj = jj+1                 !next position in long vector
             END DO
           END DO
           j = j + 1              !next slave DOF
           npsdf(j) = jj          !pointer of next slave DOF to long vector
         END DO
         !  rotates the original system onto the correct one
         IF(ndime == 2) THEN
           IF(d(1) <= 0) euler(1,ns) =  ACOS(d(2))
           IF(d(1) > 0 ) euler(1,ns) = -ACOS(d(2))
           j0 = j - 1             !slave DOF associated to rotation
           jj = npsdf(j0)         !pointer to long vector
           DO n=1,2
             alpha = COS(euler(1,ns))*l
             DO i=1,ndime
               ftsdf(jj) = alpha
               alpha      = SIN(euler(1,ns))*l
               jj = jj+1
             END DO
             l=-l
           END DO
         ELSE           ! ndime = 3
           CALL vecpro(euler(7:9,ns),d,phi)
           d = MATMUL(phi,RESHAPE(euler(1:9,ns),(/3,3/)))
           CALL actrot(euler(1:9,ns),d(1:3))
           ! slave DOF associated to rotation is j-2
           jj = npsdf(j-2)       !pointer to long vector
           kk = 3
           DO k=4,5
             DO n=1,2
               DO i=1,ndime
                 ftsdf(jj) = l*euler(kk+i,ns)
                 jj = jj+1                !next pointer to long vector
               END DO
               l  = -l
             END DO
             kk = 0
             l = -l
           END DO
         END IF     ! ndime?
       ELSE
         ifpre(ndime+1:nrotd,ns) = 1   !no DOF
       END IF     ! naeul(ns)?
     END IF     ! eul?
   END IF
 END DO     !m=1,naris
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE arisin
