      SUBROUTINE force2(isn,issdb,rssdb,lcseg,indco,cprop,x,prdat,dtime,  &
                        emass,fcont,surtf,cursl,tn,icnod,press,wear,wwear,&
                        mtdof,sldof,area,td)

!.... form contact residual force for a node-segment 2-d contact element
      USE npo_db, ONLY : tempe,tresi

      IMPLICIT NONE
!     dummy arguments
      LOGICAL, INTENT(IN) :: cursl,press,wear
      INTEGER (kind=4), INTENT (IN) :: indco,isn,lcseg(:,:),icnod,mtdof,sldof
      INTEGER (kind=4), INTENT (IN OUT) :: issdb(:)
      REAL (kind=8), INTENT (IN) :: cprop(:),x(:,:),emass(:,:), &
                     dtime,tn(:,:),td
      REAL(kind=8),POINTER :: area(:)
      REAL (kind=8), INTENT (IN OUT) :: rssdb(:),surtf(:),fcont(:,:),     &
                                        prdat(:),wwear(:)
!     local variables
      LOGICAL :: last,same
      INTEGER (kind=4) :: n1,n2,n1o,n2o,nearn,onear,n,node
      REAL (kind=8) :: fc,r,s,ro,so,to,ds,dt,pgap,pnltc,cofri,m,maux,     &
                       tslip,eslip,fricm,frict,lengt,auxil,xfact,pslip,   &
                       vns(2,3),t1(2),ts(2),eresf(2,3)

      INTERFACE
        INCLUDE 'vecuni.h'
      END INTERFACE

!.... read gap, normal direction and isoparametric coordinates of the
!.... projection point, from real-type slave surface database
      IF( prdat(1) < cprop(4) )THEN
        WRITE(58,"(' gap > cutof, node:',i6,2e15.4)",ERR=9999) isn,prdat(1),rssdb(2)
!        WRITE(* ,"(' gap > cutof, node:',i6,2e15.4)")isn,prdat(1),rssdb(2)
      END IF
      pgap     = - cprop(1) * prdat(1)*(1d0+prdat(1)/cprop(4)) !normal penalty * gap
      IF( cursl .AND. indco /= 1)THEN
        IF(indco == 0)THEN
          t1 = (-tn(:,icnod)+prdat(2:3))/2d0       !average normal
          vns(:,1) = t1/SQRT(DOT_PRODUCT(t1,t1))
        ELSE         ! IF(indco == 2)THEN
          vns(:,1) = -tn(:,icnod)              !node normal direction
        END IF
      ELSE
        vns(:,1) = prdat(2:3)          !normal direction
      END IF
                                   !segment linear shape functions
      s     = prdat(4)             !node 2 (xita)
      r     = 1d0 - s              !node 1 (1-xita)
      nearn = issdb(1)             !segment

      !  computes equivalent mass
      m = 0d0  !initializes
      maux = MAXVAL(emass(1:2,isn))   ! mass of slave node
      IF( maux /= 0d0) m = 1d0/maux
      !IF ((indco == 0) .OR. (indco == 2)) THEN  ! IF forces in both surfaces or master
      IF ( m == 0d0 ) THEN  ! IF forces in both surfaces or master
        !....   consider master nodes mass
        DO n =1,2                         !For each master node
          node = lcseg(n,nearn)           !global node number
          maux = MAXVAL(emass(1:2,node))  !max mass in any direction
          IF(maux > 0)THEN                ! add to equivalent mass
            SELECT CASE (n)
            CASE (1)                      !first node
               m = m + r**2/maux            !use r
            CASE (2)                      !second node
               m = m + s**2/maux            !use s
            END SELECT
          END IF
        END DO
      END IF
      IF(m == 0d0) THEN
        WRITE(*,    *)'  Error in nodal mass ',isn,lcseg(1:2,nearn)
        WRITE(lures,*,ERR=9999) '  Error in nodal mass, Internal numbers are: Slave ',isn, &
                              & '  Master segment ',lcseg(1:2,nearn)
        WRITE(lures,*,ERR=9999) '  Check data, at least one surface Must have mass'
        CALL runen3('CONTA2: zero masses in contact pair')
      END IF
      xfact = 2d0/(m*dtime**2)     !computes xfact
!.... form residual force vector
      IF( indco == 0 .OR. indco == 2) THEN !forces on both surfaces or master
        !....   form ns operator
        vns(:,2)   = -r*vns(:,1)    !factor for node 1
        vns(:,3)   = -s*vns(:,1)    !factor for node 2
        eresf = pgap*vns                !normal residual forces
      ELSE IF(indco == 1)THEN     !forces only in slave surface
        eresf(:,1)   = pgap*vns(:,1)  !normal residual forces
      END IF

      n1 = lcseg(1,nearn)             !nodes defining master segment
      n2 = lcseg(2,nearn)
      pslip = 0d0
!     friction treatment
      pnltc = cprop(2)                !tangential penalty parameter
      last  = issdb(2) > 0            !contact at previous step
      IF( issdb(3) < 1 ) issdb(3) = 1 !change to penetration
      IF( issdb(3) == 1)THEN          !no-slip in previous step
        cofri = cprop(3)              !Static friction coefficient
      ELSE
        cofri = cprop(6)              !Kinetic friction coefficient
      END IF
      IF(cofri /= 0D0) THEN           !if no friction, that's all
        !  computes friction forces  tslip: total slip;  eslip : elastic slip
        IF( .NOT.last ) THEN          ! IF no penetration in previous step
          rssdb(1) = s                ! stores onset natural coordinates
          frict = 0d0
        ELSE                          ! Penetration in previous step
          fricm = cofri*pgap           !maximum friction force
          t1 = x(:,n2) - x(:,n1)        !side vectors n1->n2
          so = rssdb(1)                 !old onset coordinates
          onear = issdb(2)              !previous segment
          same = onear == nearn         !same of different segments
          IF( same )THEN  ! IF penetration in previous step was in the same segment
            ds = s - so                 !differences in natural coordinates
            ts = t1*ds                  !slip vector
            CALL vecuni(2,ts,tslip)     !ts = unit vector, tslip = total slip
            frict = pnltc*tslip         !tangential force (+)
            IF(frict > fricm) THEN      !compare with maximum force
              ! IF friction force > maximum friction force ==> slip
              frict = fricm             !assign maximum force
              eslip = frict/pnltc       !maximum slip
              rssdb(1) = s - ds*eslip/tslip   !updates onset natural coordinates
              pslip = tslip - eslip     ! (+)
              issdb(3) = 2              !remember that point is sliding
            ELSE
              issdb(3) = 1              !point is fixed
            END IF
          ELSE  ! segment in previous step was in other segment
            n1o= lcseg(1,onear)     !nodes defining previous master segment
            n2o= lcseg(2,onear)
            ro = 1d0-so             !previous r coord.
            !vector between actual and onset point
            ts= r *x(:,n1) + s *x(:,n2) -ro*x(:,n1o) -so*x(:,n2o)
            CALL vecuni(2,t1,lengt)       ! side length
            auxil = DOT_PRODUCT(ts,vns(:,1))    ! projects vector over normal
            ts = ts - auxil*vns(:,1)      ! orthogonal projection
            CALL vecuni(2,ts,tslip)       ! ts = unit vector, tslip = total slip
            frict = pnltc*tslip           ! friction force
            IF(frict > fricm) THEN        ! compare with maximum friction force
              ! IF friction force > maximum friction force ==> slip
              eslip = fricm/pnltc         ! maximum relative slip
              frict = fricm               ! assign maximum friction force
              pslip = tslip - eslip
              issdb(3) = 2                !remember that point is sliding
            ELSE
              eslip = tslip               ! all tslip is elastic
              issdb(3) = 1                !point is fixed
            END IF  !maximum tangential force exceded
            auxil = DOT_PRODUCT(t1,ts)*eslip/lengt  ! slip Ds coordinate
            rssdb(1) = s - auxil        ! store onset natural coordinate
          END IF  !same or different segment
          IF( cursl .AND. indco /= 1 .AND. frict /= 0d0)THEN
            auxil = DOT_PRODUCT(ts,vns(:,1))
            ts = ts - auxil*vns(:,1)
            ts = ts/SQRT(DOT_PRODUCT(ts,ts))
          END IF
          IF( indco == 0 .OR. indco == 2) THEN !forces in both surfaces or master
            vns(1:2,1)  =  -ts(1:2)       ! slave node
            vns(1:2,2)  = r*ts(1:2)       ! master nodes
            vns(1:2,3)  = s*ts(1:2)
            eresf = eresf + frict*vns     ! add to normal forces
          ELSE IF( indco == 1 )THEN       ! forces in slave surface only
            eresf(1:2,1) = eresf(1:2,1) - frict*ts(1:2)  ! add to normal forces
          END IF
          IF ( pslip /= 0d0 ) pslip = ABS(frict*pslip*xfact)          

        END IF  !projection in previous step

        issdb(2) = nearn              !store actual segment

      END IF !if friction

      ! add to total forces if required
      surtf = surtf - eresf(1:2,1)*xfact*dtime
      IF( press ) prdat(1) = pgap*xfact*dtime

      ! Add forces into global contact forces
      IF(indco == 0 .OR. indco == 1) &             !forces on SLAVE
        fcont(1:2,isn) = fcont(1:2,isn) - eresf(1:2,1)*xfact
      IF(indco == 0 .OR. indco == 2)THEN           !forces on MASTER
        DO n=1,2                 !for each node
          node = lcseg(n,nearn)  !global number
          fcont(1:2,node) = fcont(1:2,node) - eresf(1:2,1+n)*xfact
        END DO
      END IF
      IF( wear .AND. issdb(3) == 2 )THEN
        wwear(isn) = wwear(isn) + pslip
        node = lcseg(1,nearn)  !global number
        wwear(node) = wwear(node) + pslip*r
        node = lcseg(2,nearn)  !global number
        wwear(node) = wwear(node) + pslip*s
      END IF

      !  thermal part
      IF (therm) THEN
        ! USE npo_db, ONLY : iftmp,tempe,tresi
        ! USE ctrl_db, ONLY : therm
        n = iftmp(sldof,isn)                   !  thermal dof on slave
        n1o = iftmp(mtdof,n1)                  !  thermal dof on master
        n2o = iftmp(mtdof,n2)                  !  thermal dof on master
        ! CONDUCTION        
        IF( cprop(7) > 0d0 )THEN
          fc = pgap * area(icnod)  !  force/associated area
          to = cprop(7)*(fc/cprop(9))**cprop(8)  !  convection coefficient
          so = tempe(sldof,isn)                  !  temperature on slave surface
          ro = tempe(mtdof,n1)*r + tempe(mtdof,n2)*s !master
          dt = (so - ro)*to                      !  thermal gap * h
          IF( n > 0 )   tresi(n) = tresi(n) + dt
          IF(n1o > 0) tresi(n1o) = tresi(n1o) - r*dt
          IF(n2o > 0) tresi(n2o) = tresi(n2o) - s*dt
        END IF
        ! FRICTION WORK
        IF( issdb(3) == 2 )THEN
          pslip = pslip/td !friction power 
          IF( n > 0 )  tresi(n) = tresi(n) - pslip*cprop(10)
          pslip = pslip*(1d0-cprop(10))
          IF(n1o > 0) tresi(n1o) = tresi(n1o) - pslip*r
          IF(n2o > 0) tresi(n2o) = tresi(n2o) - pslip*s
        END IF
      END IF

      RETURN
      9999 CALL runen2('')
      END SUBROUTINE force2
