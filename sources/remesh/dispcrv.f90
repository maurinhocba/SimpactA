SUBROUTINE dispcrv(npoin,nelem,lnods,lside,ldv2sd,nodcre,nwlnd,lindi,cd,cvdisp)
!-----------------------------------------------------------------------
! compute displacements of the new nodes by the curvature of the sheet
!-----------------------------------------------------------------------
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: npoin,          & !Number of nodes
                               nelem             !Number of elements in the set
  INTEGER(kind=4),POINTER:: lnods(:,:),     & !Conectivities of the set
                            ldv2sd(:,:),    & !Element side subdivision
                            lside(:,:),     & !Neighbor elements of the set
                            nodcre(:,:),    & !Nodes between is generated a new node
                            nwlnd(:,:),     & !new conectivities of refined mesh
                            lindi(:)          !label of new element by old elements
  REAL(kind=8),INTENT(IN):: cd(:,:)     !Nodal coordinates
  REAL(kind=8),POINTER:: cvdisp(:,:)    !Displacements of culvature correction

  !--- Local variables
  INTEGER(kind=4),PARAMETER:: nxn(3) = (/ 2,3,1 /) !next node in connectv.
  INTEGER(kind=4):: i, j, k, nn(3), ielem, jelem, kelem, ipoin, elcre(0:6), newnds, newel

  REAL(kind=8):: ls(3,2), t(3), n(3), long, curv, area, zdisp
  LOGICAL:: found
  INTEGER(kind=4),ALLOCATABLE:: bnod1(:), bnod2(:)
  REAL(kind=8),ALLOCATABLE:: tn(:,:), tw(:,:)
  LOGICAL,ALLOCATABLE:: actnd(:)

  ! ****  Compute nodal normals and elements Z-increments of a surface  ***
  newnds = nodcre(1,0)   !Number of new nodes in the refinement process
  newel = nodcre(2,0)    !Number of elements in the new refined mesh

  ALLOCATE(tn(3,npoin),tw(3,newnds),actnd(npoin))
  actnd(1:npoin) = .TRUE.
  tn(1:3,1:npoin) = 0d0         !nodal average normals initialization
  tw(1:3,1:newnds) = 0d0        !new nodal normals initialization

  ALLOCATE(bnod1(npoin),bnod2(npoin))   !Identify boundary nodes
  bnod1(1:npoin) = 0            !boundary  first node initialization (0 for inner node)
  bnod2(1:npoin) = 0            !boundary second node initialization (0 for inner node)

  ! Compute average normals at the nodes and the new nodal points
  kelem = 1
  found = .FALSE.
  DO ielem=1,nelem                         !for each element

    actnd(lnods(1:3,ielem)) = .FALSE.
    nn(1:3) = lnods(1:3,ielem)                   !connectivities of the element
    DO i=1,2                                     !first two side elements
      k = nxn(i)                                 !next node
      ls(1:3,i) = cd(1:3,nn(k)) - cd(1:3,nn(i))  !tang. vector
    END DO
    CALL vecpro(ls(1:3,1),ls(1:3,2),t)              !normal*area*2
    area = 0.5d0*DSQRT(DOT_PRODUCT(t(1:3),t(1:3)))  !Element area
    t(1:3) = t(1:3)/(2d0*area*area)       !Unitary normal divided by the area
    DO i=1,3
      k = nn(i)              !node number
      tn(1:3,k) = tn(1:3,k) + t  !sums on average unitary normal weighted by the inverse of the area
      IF (lside(i,ielem) /= 0) CYCLE   !Identify nodes of the boundary
      j = nn(nxn(i))    !Next node
      bnod1(k) = j      !Identify first node of the boundary
      bnod2(j) = k      !Identify second node of the boundary
    END DO

    !Found elements that created new nodal points
    IF (.NOT.found) THEN
      jelem = kelem - 1
      elcre = 0
      DO
        jelem = jelem + 1
        IF (lindi(jelem) /= ielem) CYCLE
        elcre(1) = jelem
        k = 1
        DO
          jelem = jelem + 1
          IF (jelem > newel) EXIT
          IF (lindi(jelem) == ielem) THEN
            k = k + 1                      !Se trata de un elemento divido en k veces
            elcre(k) = jelem               !Numero del k-esimo elemento de la division
          ELSE
            EXIT
          END IF
        END DO
        elcre(0) = k                   !Elemento dividido en k-veces
        IF (k > 1) found=.TRUE.
        EXIT
      END DO
      kelem = kelem + elcre(0)
    END IF

    !Compute average normals at the new nodal points
    IF (found) THEN
      found = .FALSE.
      SELECT CASE (elcre(0))
      CASE (2)
        k = nwlnd(1,elcre(2)) - npoin   !Number of node created
        tw(1:3,k) = tw(1:3,k) + t(1:3)        !sums on average unitary normal weighted by the inverse of the area
      CASE (4)
        IF (ldv2sd(3,ielem) == 0) THEN  !is a element divided by 2
          IF (lnods(2,ielem) == nwlnd(2,elcre(3))) THEN
            found = .TRUE.   !First of the two primitive elements
            k = nwlnd(1,elcre(3)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          ELSE   ! IF ( lnods(2,ielem) == nwlnd(3,elcre(4)) ) THEN
            k = nwlnd(1,elcre(4)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          END IF
        ELSE                         !is a single primitive element
          DO i=2,4
            k = nwlnd(1,elcre(i)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          END DO
        END IF
      CASE (5)
        IF ( lnods(2,ielem) == nwlnd(2,elcre(3)) ) THEN  !This is the first element if second is refined twice
          found = .TRUE.   !First of the two primitive elements
          k = nwlnd(1,elcre(3)) - npoin
          tw(1:3,k) = tw(1:3,k) + t(1:3)
        ELSE IF ( lnods(2,ielem) == nwlnd(2,elcre(5)) ) THEN  !This is the second element if first is refined twice
          DO i=4,5
            k = nwlnd(1,elcre(i)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          END DO
        ELSE IF ( lnods(2,ielem) == nwlnd(2,elcre(4)) ) THEN  !This is the first element if it is refined twice
          found = .TRUE.   !First of the two primitive elements
          DO i=4,5
            k = nwlnd(1,elcre(i)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          END DO
        ELSE ! IF ( lnods(2,ielem) == nwlnd(3,elcre(3)) ) THEN  !This is the second element if first is refined twice
          k = nwlnd(1,elcre(3)) - npoin
          tw(1:3,k) = tw(1:3,k) + t(1:3)
        END IF
      CASE (6)
        IF ( lnods(2,ielem) == nwlnd(2,elcre(3)) ) THEN  !This is the first element (refined twice)
          found = .TRUE.   !First of the two primitive elements
          DO i=3,4
            k = nwlnd(1,elcre(i)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          END DO
        ELSE ! IF ( lnods(2,ielem) == nwlnd(2,elcre(6)) ) THEN  !This is the second element (refined twice too)
          DO i=5,6
            k = nwlnd(1,elcre(i)) - npoin
            tw(1:3,k) = tw(1:3,k) + t(1:3)
          END DO
        END IF
      END SELECT
    END IF

  END DO

  ! Compute unitary normals at the nodes
  DO ipoin=1,npoin
    IF (actnd(ipoin)) CYCLE
    t(1:3) = tn(1:3,ipoin)                    !weighted normal
    long = DSQRT(DOT_PRODUCT(t(1:3),t(1:3)))  !length
    IF (long > 0d0) tn(1:3,ipoin) = t(1:3)/long  !unit vector for nodes
  END DO

  ! Compute curvatures at the new nodes
  DO i=1,newnds
    t(1:3) = tw(1:3,i)                        !weighted normal
    long = DSQRT(DOT_PRODUCT(t(1:3),t(1:3)))  !length
    IF (long > 0d0) tw(1:3,i) = t(1:3)/long   !unit vector (for new nodes)

    nn(1:2) = nodcre(1:2,i)
    found = bnod1(nn(1))==nn(2) .OR. bnod2(nn(1))==nn(2)  !TRUE for a boundary new node
    IF ( found ) THEN  !new node is a boundary node
      cvdisp(1:3,i) = 0d0
    ELSE                                         !new node is a inner node
      n(1:3) = tn(1:3,nn(2)) - tn(1:3,nn(1))  !Normal increment
      t(1:3) = cd(1:3,nn(2)) - cd(1:3,nn(1))  !Tangent vector
      long = DSQRT(DOT_PRODUCT(t(1:3),t(1:3)))    !Nodal distant
      t = t/long                                  !Unitary tangent vector
      curv = DOT_PRODUCT(n(1:3),t(1:3))           !Curvature at the new node * long

      zdisp = -0.25d0*curv*curv + 1d0
      IF (zdisp <= 0d0) zdisp=0d0
      zdisp = DSQRT(zdisp)
      zdisp = DSQRT((1d0-zdisp)/(1d0+zdisp))*long*0.5d0  !Nodal value of the displacement by curvature
      IF (curv < 0d0) zdisp=-zdisp

      cvdisp(1:3,i) = zdisp*tw(1:3,i)   !Nodal displacement by curvature
    END IF
  END DO

  DEALLOCATE(tn,tw,bnod1,bnod2,actnd)  !release memory

RETURN
END SUBROUTINE dispcrv
