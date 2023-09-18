 SUBROUTINE force3(isn,issdb,rssdb,lcseg,indco,cprop,x,prdat,dtime,  &
                   emass,eresf,cursl,tn,icnod,press,wear,lwear,      &
                   mtdof,sldof,area,td)
 !.... form contact residual force for a node-triangle 3-d contact element
 !.... thermal parte NOT PARALLELIZED YER
 USE npo_db, ONLY : tempe,tresi
 USE gvar_db, ONLY : updiv
 !$ USE omp_lib
 IMPLICIT NONE

   !--- Dummy arguments
   LOGICAL,INTENT(IN):: cursl,press,wear
   INTEGER(kind=4),INTENT(IN):: indco,isn,lcseg(:,:),icnod,mtdof,sldof
   INTEGER(kind=4),INTENT(IN OUT):: issdb(3)
   REAL(kind=8),INTENT(IN):: cprop(:),x(:,:),emass(:,:), dtime, tn(:,:),td
   REAL(kind=8),POINTER :: area(:)
   REAL(kind=8),INTENT(IN OUT):: rssdb(:),prdat(:),eresf(1:3,0:3),lwear(0:3)
   !--- Local variables
   LOGICAL:: last,same
   INTEGER(kind=4):: n1,n2,n3,n1o,n2o,n3o,nearn,onear,n,node
   REAL(kind=8):: fc, r,s,t,ro,so,to,ds,dt,pgap,pnltc,cofri,m,maux, &
                  tslip,eslip,fricm,frict,a2,auxil,xfact,pslip,     &
                  vns(3,0:3),t1(3),t2(3),t3(3),ts(3),ta(3)

   INTERFACE
     INCLUDE 'vecuni.h'
     INCLUDE 'vecpro.h'
   END INTERFACE

 !.... read gap, normal direction and isoparametric coordinates of the
 !.... projection point, from real-type slave surface database
 !  IF( prdat(1) < cprop(4) )THEN
 !    WRITE(58,"(' gap > cutof, node:',i6,2e15.4)",ERR=9999)isn,prdat(1),rssdb(3)
 !  END IF

   nearn = issdb(1)             !segment
                                !triangular area coordinates
   s     = prdat(5)             !node 2 (xita)
   t     = prdat(6)             !node 3 (eta)
   r     = 1d0 - s - t          !node 1 (zeta)
   !Computes equivalent mass
   m = 0d0  !initializes
   maux = MAXVAL(emass(1:3,isn))   ! mass of slave node
   IF (maux /= 0d0) m=1d0/maux
   !IF ((indco == 0) .OR. (indco == 2)) THEN  ! IF forces in both surfaces or master
   IF ( m == 0d0 ) THEN  ! IF forces in both surfaces or master
     !....   consider master nodes mass
     DO n=1,3                          !For each master node
       node = lcseg(n,nearn)           !global node number
       maux = MAXVAL(emass(1:3,node))  !max mass in any direction
       IF (maux > 0) THEN              ! add to equivalent mass
         SELECT CASE (n)
         CASE (1)                      !first node
           m = m + r**2/maux            !use r
         CASE (2)                      !second node
           m = m + s**2/maux            !use s
         CASE (3)                      !third node
           m = m + t**2/maux            !use t
         END SELECT
       END IF
     END DO
   END IF
   IF (m == 0d0) THEN
     WRITE(*,    *)'  Error in nodal mass ',isn,lcseg(1:3,nearn)
     WRITE(lures,*,ERR=9999)'  Error in nodal mass, Internal numbers are: Slave ',isn, &
                   '  Master segment ',lcseg(1:3,nearn)
     WRITE(lures,*,ERR=9999)'  Check data, at least one surface Must have mass'
     CALL runend('CONTA3: zero masses in contact pair')
   END IF
   xfact = 2d0/(m*dtime**2)     !computes xfact

   pgap     = - cprop(1) * prdat(1)*(1d0+prdat(1)/cprop(4)) !normal penalty * gap
   IF (cursl .AND. (indco /= 1)) THEN
     IF (indco == 0) THEN
       t1 = (-tn(:,icnod)+prdat(2:4))/2d0       !average normal
       vns(:,0) = t1/SQRT(DOT_PRODUCT(t1,t1))
     ELSE         ! IF(indco == 2)THEN
       vns(1:3,0) = -tn(1:3,icnod)              !node normal direction
     END IF
   ELSE
     vns(1:3,0) = prdat(2:4)          !normal direction
   END IF

   !.... form residual force vector
   IF ((indco == 0) .OR. (indco == 2)) THEN   !forces on both surfaces or master
     !....   form ns operator
     vns(1:3,1)   = -r*vns(1:3,0)    !factor for node 1
     vns(1:3,2)   = -s*vns(1:3,0)    !factor for node 2
     vns(1:3,3)   = -t*vns(1:3,0)    !factor for node 3
     eresf = pgap * vns     !normal residual forces
   ELSE IF(indco == 1)THEN     !forces only in slave surface
     eresf(1:3,0) = pgap * vns(1:3,0)  !normal residual forces
   END IF

   n1 = lcseg(1,nearn)             !nodes defining master segment
   n2 = lcseg(2,nearn)
   n3 = lcseg(3,nearn)
   pslip = 0d0
   !Friction treatment
   pnltc = cprop(2)                !tangential penalty parameter
   last  = issdb(2) > 0            !contact at previous step
   IF (issdb(3) < 1 .AND. updiv ) issdb(3)=1    !change to penetration
   IF (issdb(3) == 1) THEN         !no-slip in previous step
     cofri = cprop(3)              !Static friction coefficient
   ELSE
     cofri = cprop(6)              !Kinetic friction coefficient
   END IF
   IF(cofri /= 0D0) THEN           !if no friction, that's all
     !Computes friction forces  tslip: total slip;  eslip : elastic slip
     IF (.NOT.last) THEN           ! IF no penetration in previous step
       IF( updiv )THEN
         rssdb(1) = s                !stores onset natural coordinates
         rssdb(2) = t
       END IF
       frict = 0d0
     ELSE                          ! Penetration in previous step
       fricm = cofri*pgap              !maximum friction force,
       t1 = x(1:3,n2) - x(1:3,n1)    !side vectors n1->n2
       t2 = x(1:3,n3) - x(1:3,n1)    !             n1->n3
       so = rssdb(1)                 !old onset coordinates
       to = rssdb(2)
       onear = issdb(2)              !previous segment
       same = onear == nearn         !same of different segments
       IF (same) THEN  ! IF penetration in previous step was in the same segment
         ds = s - so                 !differences in natural coordinates
         dt = t - to
         ts = t1*ds + t2*dt          !slip vector
         CALL vecuni(3,ts,tslip)     !ts = unit vector, tslip = total slip
         frict = pnltc*tslip         !tangential force
         IF (frict > fricm) THEN     !compare with maximum force
           ! IF friction force > maximum friction force ==> slip
           frict = fricm             !assign maximum force
           eslip = frict/pnltc       !maximum slip
           IF( updiv )THEN
             rssdb(1) = s - ds*eslip/tslip   !updates onset natural coordinates
             rssdb(2) = t - dt*eslip/tslip
           END IF
           pslip = tslip - eslip
           IF( updiv )issdb(3) = 2              !remember that point is sliding
         ELSE
           IF( updiv )issdb(3) = 1              !point is fixed
         END IF
       ELSE  ! segment in previous step was in other segment
         n1o = lcseg(1,onear)    !nodes defining previous master segment
         n2o = lcseg(2,onear)
         n3o = lcseg(3,onear)
         ro = 1d0-so-to          !previous r coord.
         !vector between actual and onset point
         ts= r *x(:,n1) + s *x(:,n2) + t *x(:,n3)   &
            -ro*x(:,n1o)- so*x(:,n2o)- to*x(:,n3o)
         CALL vecpro(t1,t2,t3)         ! t3 = t1 x t2
         CALL vecuni(3,t3,a2)          ! a2 = twice the area
         auxil = DOT_PRODUCT(ts,t3)    ! projects vector over normal
         ts = ts - auxil*t3            ! orthogonal projection
         CALL vecuni(3,ts,tslip)       ! ts = unit vector, tslip = total slip
         frict = pnltc*tslip           ! friction force
         IF(frict > fricm) THEN        ! compare with maximum friction force
           ! IF friction force > maximum friction force ==> slip
           eslip = fricm/pnltc         ! maximum relative slip
           frict = fricm               ! assign maximum friction force
           pslip = tslip - eslip
           IF( updiv )issdb(3) = 2                !remember that point is sliding
         ELSE
           eslip = tslip               ! all tslip is elastic
           IF( updiv )issdb(3) = 1                !point is fixed
         END IF  !maximum tangential force exceded
         CALL vecpro(ts,t2,ta)       ! ta = ts x t2
         auxil = DOT_PRODUCT(ta,t3)/a2 *eslip  ! slip Ds coordinate
         IF( updiv )rssdb(1) = s - auxil        ! store onset natural coordinate
         CALL vecpro(t1,ts,ta)       ! ta = t1 x ts
         auxil = DOT_PRODUCT(ta,t3)/a2 *eslip  ! slip Dt coordinate
         IF( updiv )rssdb(2) = t - auxil        ! store onset natural coordinate
       END IF  !same or different segment
       IF (cursl .AND. (indco /= 1) .AND. (frict /= 0d0)) THEN
         a2 = DOT_PRODUCT(ts,vns(:,0))          !orthogonal projection
         ts = ts - a2*vns(:,0)
         ts = ts/SQRT(DOT_PRODUCT(ts,ts))
       END IF
       IF ((indco == 0) .OR. (indco == 2)) THEN  !forces in both surfaces or master
         vns(1:3,0)  =  -ts(1:3)       ! slave node
         vns(1:3,1)  = r*ts(1:3)       ! master nodes
         vns(1:3,2)  = s*ts(1:3)
         vns(1:3,3)  = t*ts(1:3)
         eresf = eresf + frict*vns     ! add to normal forces
       ELSE IF( indco == 1 )THEN       ! forces in slave surface only
         eresf(1:3,0) = eresf(1:3,0) - frict*ts(1:3)  ! add to normal forces
       END IF
       IF ( pslip /= 0d0 ) pslip = ABS(frict*pslip*xfact)

     END IF  !projection in previous step

     IF( updiv )issdb(2) = nearn              !store actual segment

   END IF !if friction

   ! Add to total forces if required
   IF (press) prdat(1)=pgap*xfact*dtime
   SELECT CASE (indco)
   CASE (0)               !forces on BOTH
     eresf =  + eresf*xfact
   CASE (1)               !forces on SLAVE
    eresf(1:3,0) =  + eresf(1:3,0)*xfact
   CASE (2)               !forces on MASTER
     eresf(1:3,1:3) = + eresf(1:3,1:3)*xfact
   END SELECT

   ! add to friction work (wearing measure)
   IF (wear .AND. issdb(3) == 2 .AND. updiv )  lwear = + pslip*(/ 1d0,r,s,t /)

   !  thermal part
   IF (therm) THEN
     !! USE npo_db, ONLY : iftmp,tempe,tresi
     !! USE ctrl_db, ONLY : therm
     ! CONDUCTION
     n = iftmp(sldof,isn)                   !  thermal dof on slave
     n1o = iftmp(mtdof,n1)                  !  thermal dof on master
     n2o = iftmp(mtdof,n2)                  !  thermal dof on master
     n3o = iftmp(mtdof,n3)                  !  thermal dof on master
     IF( cprop(7) > 0d0 )THEN
       fc = pgap * area(icnod)  !  force/associated area
       to = cprop(7)*(fc/cprop(9))**cprop(8)  !  convection coefficient
       so = tempe(sldof,isn)                  !  temperature on slave surface
       ro = tempe(mtdof,n1)*r + tempe(mtdof,n2)*s + tempe(mtdof,n3)*t !master
       dt = (so - ro)*to                      !  thermal gap * h
       IF( n > 0 )  tresi(n) = tresi(n) + dt
       IF(n1o > 0) tresi(n1o) = tresi(n1o) - r*dt    !lock TRESI (NOT YET)
       IF(n2o > 0) tresi(n2o) = tresi(n2o) - s*dt
       IF(n3o > 0) tresi(n3o) = tresi(n3o) - t*dt
     END IF
     !  FRICTION WORK
     IF( issdb(3) == 2 )THEN
       pslip = pslip/td    !friction work per unit time
       IF( n > 0 )  tresi(n) = tresi(n) - pslip*cprop(10)    !lock TRESI
       pslip = pslip*(1d0-cprop(10))
       IF(n1o > 0) tresi(n1o) = tresi(n1o) - pslip*r
       IF(n2o > 0) tresi(n2o) = tresi(n2o) - pslip*s
       IF(n3o > 0) tresi(n3o) = tresi(n3o) - pslip*t
     END IF
   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE force3
