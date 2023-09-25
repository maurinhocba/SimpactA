 SUBROUTINE stra11(ntype,ngaus,x,ni,nd,fix,ln, &
                   l0,cab,r0,lambd,pff,curv,t,n,sh)
 !***********************************************************************
 !
 !*****this routine computes strain measures for 2-d beam/shell elements
 !
 !***********************************************************************
 IMPLICIT NONE
 !              routine parameters
 INTEGER (kind=4), INTENT(IN) :: ntype,  & !problem type
                                 ngaus,  & !number of integration points
                                 ln(4)     !1-4 indicate branching nodes or boundary
 LOGICAL, INTENT(IN) :: fix(:)             !(2) TRUE: constrained  FALSE:free
 REAL (kind=8), INTENT(IN) :: x(2,4),  & ! node coordinates
                              sh(2,2), & ! shape functions
                              l0(3),   & ! initial length
                              r0(2,2)      ! initial X1 coordinate at Gauss points
 REAL (kind=8), INTENT(IN OUT) :: cab(4,2)   ! coefficients for curvature evaluation and others
 REAL (kind=8), INTENT(OUT) :: pff(:,:),      & ! stretchings (squared)
                               curv(:,:),     & ! curvatures (twice)
                               lambd(3),      & ! stretchings
                               t(2),          & ! tangent vector
                               n(2)             ! normal vector


 TYPE (nodei), POINTER :: ni,nd           ! pointers to branching nodes

 !              local variables
 INTEGER (kind=4) :: g
 REAL (kind=8) r,li,h(2),ro,gamma,c1,c2,ca,sa

 ! basic computations
 t   = x(:,3) - x(:,2)               !tangent vector
 li  = SQRT(DOT_PRODUCT(t,t))        !actual length
 t   = t/li                          !unit tangent vector
 n   = (/ -t(2), t(1) /)             !unit normal vector
 lambd(1) = li/l0(2)                 !stretching (s)

 !  compute angle change (gamma) at first node
 SELECT CASE (ln(1)) !left node
 CASE (:-1)          !branching node
   gamma = 2d0*(ni%ba-ATAN3(t(2),t(1),cab(1,1)))  !twice the rotation diff
   cab(2,1) = 2d0*lambd(1)/l0(2)                  !twice the inverse length (factor)
 CASE (0)            !end node
   IF( fix(1) )THEN  !clamped or symmetry plane
     gamma = -2d0*ATAN3(t(2),t(1),cab(1,1))       !twice the bar rotation
     cab(2,1) = 2d0*lambd(1)/l0(2)                !factor
   ELSE              !free or simple supported
     gamma = 0d0     !0
   END IF
 CASE (1:)           !inner node
   h = x(:,2) - x(:,1)          !tangent at side
   ro = SQRT(DOT_PRODUCT(h,h))  !length
   h = h/ro                     !unit tangent at side
   ca = DOT_PRODUCT(t,h)        !cos a
   sa = DOT_PRODUCT(n,h)        !sin a
   gamma = cab(4,1)*ATAN3(sa,ca,cab(1,1))     !distorion
   cab(2,1) = 1d0/ro            !factor (1/Li)
 END SELECT
 cab(3,1) = gamma               !keep distorsion
 c1 = lambd(1)/l0(2)*gamma               !curvature Lamba_s/L0 * G

 !  compute angle change (gamma) at second node
 SELECT CASE (ln(4)) !left node
 CASE (:-1)          !branching node
   gamma = 2d0*(ATAN3(t(2),t(1),cab(1,2))-nd%ba)  !twice the rotation diff
   cab(2,2) = 2d0*lambd(1)/l0(2)                  !twice the inverse length (factor)
 CASE (0)            !end node
   IF( fix(2) )THEN  !clamped or symmetry plane
     gamma = 2d0*ATAN3(t(2),t(1),cab(1,2))        !twice the bar rotation
     cab(2,2) = 2d0*lambd(1)/l0(2)                !factor
   ELSE              !free or simple supported
     gamma = 0d0     !0
   END IF
 CASE (1:)           !inner node
   h = x(:,4) - x(:,3)          !tangent at side
   ro = SQRT(DOT_PRODUCT(h,h))  !length
   h = h/ro                     !unit tangent at side
   ca = DOT_PRODUCT(t,h)        !cos a
   sa = -DOT_PRODUCT(n,h)       !sin a
   gamma = cab(4,2)*ATAN3(sa,ca,cab(1,2))  !distorion
   cab(2,2) = 1d0/ro            !factor
 END SELECT
 cab(3,2) = gamma             !keep distorsion
 c2  = lambd(1)/l0(2)*gamma            !curvature Lamba_s/L0 * G

 ! compute strains at Gauss points

 li = 0d0             !initializes Transverse stretching
 DO g=1,ngaus         !for each Gauss Point
   !       strains
   pff(1,g) = lambd(1)*lambd(1)              !Css (zeta = 0)
   curv(1,g) = 2d0*(sh(1,g)*c1 + sh(2,g)*c2) !twice the in-plane curvature
   SELECT CASE (ntype)
   CASE (1)    !beam
     li         = li + 1d0/SQRT(lambd(1))    !uniform
   CASE (2)    !plane strain
     pff(2,g)   = 1d0                        !1d0
     curv(2,g)  = 0d0                        !no curvature
     li         = li + 1d0                   !1d0
   CASE (3)    !axilsymmetric
     r    = sh(1,g)*x(1,2) + sh(2,g)*x(1,3)  !present x1 coordinate
     pff(2,g) = (r/r0(1,g))**2               !C33 (zeta = 0)
     curv(2,g)= - 2d0*(t(2)-r0(2,g))*r/r0(1,g)**2      !twice the second curvature
     li = li + r/r0(1,g)                     !transverse stretching
   END SELECT
 END DO
 IF( ngaus == 1 ) curv(1,2) = c2 - c1        !stabilization curvature
 lambd(3) = li/ngaus                 !transverse stretching
 lambd(2) = 1d0/lambd(1)/lambd(3)    !normal stretching
 RETURN
 END SUBROUTINE stra11
