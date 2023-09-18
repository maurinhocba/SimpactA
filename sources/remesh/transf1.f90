SUBROUTINE transf1(ne_new,nnod,ng,coorn,coora,l_old,l_new,sg_new,strnd)

  ! transfers data from the old mesh to the new one (from old nodes to new Gauss points)

  USE lispa0
  USE ctrl_db, ONLY : ndime
  IMPLICIT NONE

  ! dummy arguments
  INTEGER (kind=4), INTENT(IN) :: ne_new, &       ! number of elements in the new mesh
                                  nnod,   &       ! number of nodes per elemen (new)
                                  ng              ! number of Gauss points
  INTEGER (kind=4), POINTER    :: l_old(:,:),  &  ! connectivities of old elements
                                  l_new(:,:)      ! connectivities of new elements
  REAL (kind=8), INTENT(IN)    :: coorn(:,:),  &  ! coordinates of new points
                                  coora(:,:)      ! coordinates of old points

  REAL (kind=8), POINTER :: &
                      sg_new(:,:),  &  ! internal variables at gauss points (new mesh)
                      strnd(:,:)       ! internal variables at nodes (old mesh)

  ! local
  INTEGER (kind=4) :: i, j, k, ie, ip, lj, ne_old, nepe, nnode, &
                      nvarg, l(3,2),io,np,ig,ng1
  REAL (kind=8) :: xn(ndime), xx(ndime,nnod),x(ndime,3), ar(3), area, det2, tol, tiny
  REAL (kind=8) :: posgp(2), shape(nnod,ng),deriv(nnod,2),xi,eta
  REAL (kind=8), PARAMETER ::  SMALL=5.0D-2    !tolerance in proyection
  REAL (kind=8), PARAMETER ::  tol0 =1.0D-4    !tolerance in bounding-box
  REAL (kind=8), POINTER ::  ftot(:)

  INTERFACE
    INCLUDE 'nonear.h'
  END INTERFACE

  nnode  = SIZE(l_old,1)       !number of nodes per element in previous mesh
  ne_old = SIZE(l_old,2)       !number of elements in previous mesh
  nvarg  = SIZE(strnd,1)       !number of internal variables to transfer
  ALLOCATE (ftot(nvarg))
  IF (nnod == 4) THEN
    posgp = (/ -0.577350269189626D+00, 0.577350269189626D+00 /)
    ng1 = 2
  ELSE !IF (nnod == 3) THEN
    posgp = (/  0.333333333333333D+00, 0d0  /)
    ng1 = 1
  END IF

  !      gauss points shape and derivatives of nodal functions

  k = 0
  DO i=1,ng1
    xi = posgp(i)
    DO j=1,ng1
      k = k+1
      eta = posgp(j)
      CALL shape3(deriv(1,1),shape(1,k),xi,eta,nnod)
    END DO
  END DO

  !       reserve space for Gauss variables, velocities and auxiliar array
  ! initialzes tolerances in bounding box and local area coordinates
  tol = MAXVAL(coorn(1,:)) - MINVAL(coorn(1,:))
  tol = MAX(MAXVAL(coorn(2,:)) - MINVAL(coorn(2,:)),tol)*tol0 !bounding box tolerance
  tiny = small                                                !local coord tolerance
  ! loop over new elements

  ie = 1         !initializes first element to look in
  ig = 0         !initializes New Gauss point pointer
  DO ip = 1, ne_new       !for each new element
    xx = coorn(:,l_new(:,ip))  !coordinates
    L10: DO k=1,ng                    !for each New gauss point
      xn = MATMUL(xx,shape(:,k)) ! coordinates of new gauss points
      ig = ig+1                  ! gauss point order
      ! loop over all the elements in the old mesh to find its proyection
      io = ie                  !remember first searched element
15    np = 0                   !number of elements with proyection
      L15: DO                  !loop until proyection found

        CALL trian(nnode,l_old(1:nnode,ie),l,nepe)   !divide Quadrilateral element in Triangle

        L20: DO i=1,nepe  ! number of sub-triangles in an element (1 for Triangles)

          x(:,:) = coora(:,l(:,i))                   !triangle coordinates

          ! discard an element if this point is outside its boundingbox
          IF (nonear (xn, x, ndime,tol)) CYCLE L20
          np = np+1
          ! element area
          area  = det2(x(1,1),x(1,2),x(1,3),ndime,.TRUE.)

          ! three subelements (natural triangular coordinates)
          ar(1) = det2(x(1,2), x(1,3), xn, ndime,.FALSE.) / area
          IF (ar(1) < -TINY) CYCLE L20

          ar(2) = det2(x(1,3), x(1,1), xn, ndime,.FALSE.) / area
          IF (ar(2) < -TINY) CYCLE L20

          ar(3) = 1d0 -ar(1) -ar(2)
          !ar(3) = det2(x(1,1), x(1,2), xn, ndime,.FALSE.) / area
          IF (ar(3) < -TINY) CYCLE L20

          ! interpolate using area coordinates
          ftot = 0d0           !initializes variables
          DO j=1,3             !for each node
            lj = l(j,i)        !associated node
            ftot(:) = ftot(:) + strnd(:,lj)*ar(j)  !internal variables
          END DO
          sg_new(:,ig) = ftot(:)
          CYCLE L10  ! jump out - don't check any more elements

        END DO L20   !loop over individual triangles
        ie = MOD(ie,ne_old) + 1  !new element (if passed end of the list, back to first)
        IF( ie == io )EXIT       !if back to the first element tested =>EXIT
      END DO L15
      WRITE(lures,*,ERR=9999) ip,'  not placed, FATAL'   !proyection not found => ERROR
      WRITE(55,*,ERR=9999) ip,'  not placed, FATAL'      !proyection not found => ERROR
      PRINT *, ip,'  not placed, FATAL'                  !proyection not found => ERROR
      IF( np == 0 )THEN
        tol = tol*2d0  !increase bounding-box tolerance
      ELSE IF ( np == 1 ) THEN
        tol = tol*2d0  !increase bounding-box tolerance
        tiny= 2d0*tiny !increase local coordinate tolerance
      ELSE
        tiny= 2d0*tiny !increase local coordinate tolerance
      END IF
      GO TO 15   !loop again with a larger tolerance
      !STOP
    END DO L10
  END DO   ! end of new point loop
RETURN
 9999 CALL runen2('')
END SUBROUTINE transf1
