      SUBROUTINE force1(isn,issdb,rssdb,lcseg,indco,cprop,x,prdat,dtime,  &
                        emass,fcont,velnp,surtf,press)  !,cursl,tn

!.... form contact residual force for a tube-skin 3-d contact element

      IMPLICIT NONE
!     dummy arguments
      LOGICAL, INTENT(IN) :: press  !cursl
      INTEGER (kind=4), INTENT (IN) :: indco,isn,lcseg(:,:)
      INTEGER (kind=4), INTENT (IN OUT) :: issdb(:)
      REAL (kind=8), INTENT (IN) :: cprop(:),x(:,:),emass(:,:),velnp(:,:), &
                     dtime !,tn(:,:)
      REAL (kind=8), INTENT (IN OUT) :: rssdb(:),surtf(:),fcont(:,:),     &
                                        prdat(:)
!     local variables
      LOGICAL :: last
      INTEGER (kind=4) :: n1,n2,nearn,onear,n,node
      REAL (kind=8) :: r,s,so,ds,pgap,pnltc,cofri,m,maux,ratio,            &
                       eslip,fricm,frict,frics,fricr,lengt,xfact,rslip,tslip,    &
                       vns(3,3),ts(3),eresf(3,3),torque

      INTERFACE
        INCLUDE 'vecuni.h'
      END INTERFACE

!.... read gap, normal direction and isoparametric coordinates of the
!.... projection point, from real-type slave surface database
      pgap     = - cprop(1) * prdat(1)  !normal penalty * gap
      !IF( cursl .AND. indco /= 1)THEN
      !  IF(indco == 0)THEN
      !    t1 = (-tn(:,icnod)+prdat(2:3))/2d0       !average normal
      !    vns(:,1) = t1/SQRT(DOT_PRODUCT(t1,t1))
      !  ELSE         ! IF(indco == 2)THEN
      !    vns(:,1) = -tn(:,icnod)              !node normal direction
      !  END IF
      !ELSE
        vns(:,1) = prdat(3:5)      !normal direction
      !END IF
                                   !segment linear shape functions
      s     = prdat(2)             !node 2 (xita)
      r     = 1d0 - s              !node 1 (1-xita)
      nearn = issdb(1)             !segment

      !  computes equivalent mass
      m = 0d0  !initializes
      maux = MAXVAL(emass(1:3,isn))   ! mass of slave node
      IF( maux /= 0d0) m = 1d0/maux
      IF ( m == 0d0 ) THEN  ! IF forces in both surfaces or master
        !....   consider master nodes mass
        DO n =1,2                         !For each master node
          node = lcseg(n,nearn)           !global node number
          maux = MAXVAL(emass(1:3,node))  !max mass in any direction
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
        CALL runen3('CONTA1: zero masses in contact pair')
      END IF
      xfact = 2d0/(m*dtime**2)     !computes xfact
!.... form residual force vector
      IF( indco == 0 .OR. indco == 2) THEN !forces on both surfaces or master
        !....   form ns operator
        vns(:,2)   = -r*vns(:,1)    !factor for node 1
        vns(:,3)   = -s*vns(:,1)    !factor for node 2
        eresf = pgap*vns            !normal residual forces
      ELSE IF(indco == 1)THEN     !forces only in slave surface
        eresf(:,1)   = pgap*vns(:,1)  !normal residual forces
      END IF

!     friction treatment
      pnltc = cprop(2)                !tangential penalty parameter
      last  = issdb(2) > 0            !contact at previous step
      ! FRICTION DUE TO AXIAL ROTATION MISSING
      IF( issdb(3) == 0)THEN          !no-slip in previous step
        cofri = cprop(3)              !Static friction coefficient
      ELSE
        cofri = cprop(4)              !Kinetic friction coefficient
      END IF
      IF(cofri /= 0D0) THEN           !if no friction, that's all
        rslip = -velnp(4,isn)*cprop(6)/2d0*dtime
        !  computes friction forces  tslip: total slip;  eslip : elastic slip
        s  = s*prdat(9)               ! onset coordinate
        IF( .NOT.last ) THEN          ! IF no penetration in previous step
          rssdb(1) = s                ! stores onset coordinate
          frics = 0d0; fricr = 0d0    !no friction forces 
          issdb(2) = nearn
        ELSE                          ! Penetration in previous step
          fricm = cofri*pgap            !maximum friction force,
          ts = prdat(6:8)               !side vectors n1->n2
          so = rssdb(1)                 !old onset coordinates
          onear = issdb(2)              !previous segment
          IF( onear /= nearn )THEN  ! IF penetration in previous step was in the same segment
            IF( onear < nearn )THEN  !previous segment is before
              n1= lcseg(1,onear)       !nodes defining previous master segment
              n2= lcseg(2,onear)
              ts= - x(:,n1) + x(:,n2)
              lengt = SQRT(DOT_PRODUCT(ts,ts))   ! segment length
              so = -lengt + so
            ELSE                      !previous segment is after
              so = prdat(9) + so
            END IF
            issdb(2) = nearn
          END IF  !same or different segment
          ds = s - so                 !differences in coordinates
          tslip = SQRT(rslip**2 + ds**2) !(+)
          ratio = ds/tslip
          frict = pnltc*tslip         !tangential force (+)
          IF(frict > fricm) THEN      !compare with maximum force
            ! IF friction force > maximum friction force ==> slip
            frict = fricm             !assign maximum tangential force(+)
            issdb(3) = 2              !remember that point is sliding
            eslip = frict/pnltc*ratio !maximum slip (+/-)
            rssdb(1) = s - eslip      !updates onset coordinate
          ELSE
            issdb(3) = 1              !point is fixed
          END IF
          frics = frict*ratio         !assign maximum tangential force(+/-)
          fricr = frict*rslip/tslip   !assign maximum tangential force(+/-)
          !IF( cursl .AND. indco /= 1 .AND. frict /= 0d0)THEN
          !  auxil = DOT_PRODUCT(ts,vns(:,1))
          !  ts = ts - auxil*vns(:,1)
          !  ts = ts/SQRT(DOT_PRODUCT(ts,ts))
          !END IF
          IF( indco == 0 .OR. indco == 2) THEN !forces in both surfaces or master
            vns(1:3,1)  =  -ts(1:3)       ! slave node
            vns(1:3,2)  = r*ts(1:3)       ! master nodes
            vns(1:3,3)  = s*ts(1:3)
            eresf = eresf + frics*vns     ! add to normal forces
          ELSE IF( indco == 1 )THEN       ! forces in slave surface only
            eresf(1:3,1) = eresf(1:3,1) - frics*ts(1:3)  ! add to normal forces
          END IF

        END IF  !projection in previous step
      ELSE
        fricr = 0d0
      END IF !if friction

      ! add to total forces if required
      surtf = surtf - eresf(1:3,1)*xfact*dtime
      IF( press ) prdat(1) = pgap*xfact*dtime

!     Add forces into global contact forces

      IF(indco == 0 .OR. indco == 1) &             !forces on SLAVE
        fcont(1:3,isn) = fcont(1:3,isn) - eresf(1:3,1)*xfact
      IF(indco == 0 .OR. indco == 2)THEN           !forces on MASTER
        DO n=1,2                 !for each node
          node = lcseg(n,nearn)  !global number
          fcont(1:3,node) = fcont(1:3,node) - eresf(1:3,1+n)*xfact
        END DO
      END IF
      IF( fricr /= 0d0 )THEN
        torque = fricr*xfact*cprop(6)/2d0            !maximum torque
        IF(indco == 0 .OR. indco == 1) &             !forces on SLAVE
          fcont(4,isn) = fcont(4,isn) - torque
        IF(indco == 0 .OR. indco == 2)THEN           !forces on MASTER
          DO n=1,2                 !for each node
            node = lcseg(n,nearn)  !global number
            fcont(4,node) = fcont(4,node) + torque*r
            r = 1d0 - r
          END DO
        END IF
      END IF

      RETURN
      9999 CALL runen2('')
      END SUBROUTINE force1
