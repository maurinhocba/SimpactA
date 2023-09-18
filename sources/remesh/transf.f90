SUBROUTINE transf(numpn, coorn, coora, l_old, v_old, v_new, elemp)

  ! transfers data from the old mesh to the new one (from old to new nodes)

  USE lispa0
  USE ctrl_db, ONLY : ndime, ndofn
  IMPLICIT NONE

  ! dummy arguments
  INTEGER (kind=4) :: numpn,     &    ! number of new points
                      l_old(:,:)      ! connectivities of old elements
  INTEGER (kind=4), POINTER, OPTIONAL :: elemp(:) !element where projects
  REAL (kind=8)    :: coorn(:,:),  &  ! coordinates of new points
                      coora(:,:),  &  ! coordinates of old points
                      v_old(:,:)      ! nodal velocities (old mesh)
  REAL (kind=8), POINTER :: &
                      v_new(:,:)      ! nodal velocities (new mesh)

  ! local
  LOGICAL store           ! store element where node projects
  INTEGER (kind=4) :: i, j, ie, ip, lj, ne_old, nepe, nnode, &
                      l(3,2),io,np
  REAL (kind=8) :: xn(ndime), v(ndofn), x(ndime,3), ar(3), area, det2, tol, tiny
  REAL (kind=8), PARAMETER ::  SMALL=5.0D-2    !tolerance in proyection
  REAL (kind=8), PARAMETER ::  tol0 =1.0D-4    !tolerance in bounding-box

  INTERFACE
    INCLUDE 'nonear.h'
  END INTERFACE

  nnode  = SIZE(l_old,1)       !number of nodes per element in previous mesh
  ne_old = SIZE(l_old,2)       !number of elements in previous mesh
  !       reserve space for Gauss variables, velocities and auxiliar array
  ALLOCATE( v_new(ndofn,numpn))
  store =  PRESENT(elemp)                 !store projection?
  IF (store) ALLOCATE(elemp(numpn))
  v_new = 0d0    !initializes velocities
  ! initialzes tolerances in bounding box and local area coordinates
  tol = MAXVAL(coorn(1,:)) - MINVAL(coorn(1,:))
  tol = MAX(MAXVAL(coorn(2,:)) - MINVAL(coorn(2,:)),tol)*tol0 !bounding box tolerance
  tiny = small                                                !local coord tolerance
  ! loop over new points

  ie = 1                    !initializes first element to look in
  L10: DO ip = 1, numpn      !for each new point
     xn  = coorn(:,ip)        !new coordinates

    ! loop over all the elements in the old mesh to find its proyection
    io = ie                  !remember first searched element
15   np = 0                   !number of elements with proyection
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
        v = 0d0              !initializes velocities
        DO j=1,3             !for each node
          lj = l(j,i)        !associated node
          v(:) = v(:) + v_old(:,lj)*ar(j)        !velocity
        END DO
        v_new(:,ip) = v(:)
        IF (store) elemp(ip) = ie
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

  END DO L10  ! end of new point loop
RETURN
 9999 CALL runen2('')
END SUBROUTINE transf
