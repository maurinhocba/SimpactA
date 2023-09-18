SUBROUTINE plines(cdbnd,nlbnd,fside,ifx_b)
! Analyze contour fixity codes to find polylines
USE ifx_db,ONLY: nd1       !number of fixity codes per node
USE meshmo_db,ONLY: pline, headp, tailp, nstart, lstrt, nplins, &
                    new_pline, add_pline, r_elm_zone
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: ifx_b(:,:)   !(nd1,nelbnd) fixities at second nodes
  REAL(kind=8):: cdbnd(:,:)      !(ndime,nelbnd) coords at second nodes
  INTEGER(kind=4):: nlbnd(:)     !(nline) number of segments in each boundary line
  INTEGER(kind=4):: fside(:)     !(nelbnd) free-side segments
  !--- Local variables
  INTEGER(kind=4) :: i,i1,i2,i3,l,n1,n2,n3,nl,ip,nsegs,&
                     nfirst,nlast,frstn,frsts,segbnd
  TYPE (pline), POINTER :: pln

  INTERFACE
    INCLUDE 'vertex.h'
  END INTERFACE
  !------------------------------------------------------------

  NULLIFY(headp,tailp)   !Initializes empty list of boundary lines
  nstart = 1             !Initializes first position of a line
  nl = SIZE(nlbnd,1)     !Number of lines in boundary
  ALLOCATE(lstrt(nl))    !starting node number in each boundary line

  !loop over each line in boundary
  nplins = 0            !initializes: number of polylines
  frstn = 1             !initializes: lower nodes limits
  frsts = 0             !initializes: upper nodes limits

  DO l = 1,nl
    segbnd = frsts + nlbnd(l) !number of segments in line
    lstrt(l) = frstn          !first segment in contour line

    IF (.NOT.r_elm_zone) THEN
    ! in complete domain remeshing: creates a polylines (or NURBS)
    ! based on boundary conditions or angle between subsequents lines

      ! search for a starting point (first point with non-zero codes)
      i1 = segbnd               !point to previous node to present
      i2 = frstn                !point to present (first) possible node
      i3 = i2 + 1               !point to next node to present (NOT used)

      nfirst = frstn   !first node in the polyline (candidate)
      ip = frstn + 1   !for use in position respect to nstart
      nstart = frstn   !start in first node
      nsegs = 0

      DO      !loop to obtain the polylines
        nlast = 0      !Initializes last position of a line
        n1 = MOD((i1+nstart-ip),segbnd) + 1 !previous position respect to nstart
        n2 = MOD((i2+nstart-ip),segbnd) + 1  !present position respect to nstart
        n3 = MOD((i3+nstart-ip),segbnd) + 1   !next position respect to nstart

        ! compare boundary conditions of present point and next point
        ! to detect end of polyline (change of boundary conditions)

        ! since GiD performing Collapse operation looses boundary conditions,
        ! additional analysis of the angles have been added here
        ! and later Collapse will not be performed

        IF ( vertex(cdbnd, n1, n2, n3)               .OR. &  ! angle > 15   or
             ANY(ifx_b(1:nd1,n2) /= ifx_b(1:nd1,n3)) .OR. &  ! change of BC or
             n2 == segbnd ) THEN  ! last segment
          CALL new_pline(pln)   !allocate a new boundary line
          pln%nfirst = nfirst   !associate first position
          IF( vertex(cdbnd, n1, n2, n3) .AND. n1 /= segbnd )THEN !if vertex then
            nlast = n1 !initialize
            IF   (ANY(ifx_b(1:nd1,n1) /= 0) .AND. &
                  ANY(ifx_b(1:nd1,n2) /= 0))THEN
              pln%ifix(1:nd1) = 0 !initialize
              DO i = 1,nd1  !assign only if line have same restrictions in nodes
                IF(ifx_b(i,n1) == ifx_b(i,n2)) pln%ifix(i) = ifx_b(i,n1)
              END DO
            ELSE
              pln%ifix(1:nd1) = 0
            END IF
          ELSE ! if not vertex -> is a change of BC
            nlast = n2 !initialize
            IF   (ANY(ifx_b(1:nd1,n2) /= 0) .AND. &
                  ANY(ifx_b(1:nd1,n3) /= 0))THEN
              pln%ifix(1:nd1) = 0 !initialize
              DO i = 1,nd1  !assign only if line have same restrictions in nodes
                IF(ifx_b(i,n2) == ifx_b(i,n3)) pln%ifix(i) = ifx_b(i,n2)
              END DO
            ELSE
              pln%ifix(1:nd1) = 0
            END IF
          END IF
          nsegs = nsegs + (nlast - nfirst) + 1
          pln%nlast  = nlast                   !last node in the line
          nplins = nplins + 1                  !increase number of lines
          CALL add_pline(pln,headp,tailp)      !add to the list
          nfirst = nlast + 1                   !new first node for next line
          i1 = nfirst                          !new previous node
          i2 = i1 + 1                          !new present node
          i3 = i2 + 1                          !new next node
        ELSE           !same boundary line
          i1 = i2                              !new previous node
          i2 = i1 + 1                          !new present node
          i3 = i2 + 1                          !new next node
        END IF
        !if last segment is forward a vertex point
        IF(n2    == segbnd           .AND. &
           vertex(cdbnd, n1, n2, n3) .AND. &
           nsegs == segbnd - 1 )THEN
          CALL new_pline(pln)   !allocate a new boundary line
          pln%nfirst = nfirst   !associate first position
          nlast = n2 !initialize
          IF   (ANY(ifx_b(1:nd1,n2) /= 0) .AND. &
                ANY(ifx_b(1:nd1,n3) /= 0))THEN
            pln%ifix(1:nd1) = 0 !initialize
            DO i = 1,nd1  !assign only if line have same restrictions in nodes
              IF(ifx_b(i,n2) == ifx_b(i,n3)) pln%ifix(i) = ifx_b(i,n2)
            END DO
          ELSE
            pln%ifix(1:nd1) = 0
          END IF
          nsegs = nsegs + (nlast - nfirst) + 1
          pln%nlast  = nlast                   !last node in the line
          nplins = nplins + 1                  !increase number of lines
          CALL add_pline(pln,headp,tailp)      !add to the list
        END IF
        IF (nsegs >= nlbnd(l)) EXIT ! all segments was processed
      END DO

    ELSE
    ! in zone remeshing each segment is a polyline, this is for avoid problems
    ! with nodes position when paste patch of elements in original mesh

      ! Loop to read data and add them to the list
      i1 = frstn            !point previous to present node
      i2 = i1 + 1           !present point
      i3 = i2 + 1           !point next to present
      ! ...
      nstart = frstn        !start in first segment
      nfirst = nstart       !first node in the polyline (candidate)
      ip = frstn + 1        !for use in position respect to nstart

      DO      !loop to obtain the polylines
        n1 = MOD((i1+nstart-ip),segbnd) + 1   !previous position respect to nstart
        n2 = MOD((i2+nstart-ip),segbnd) + 1   !present position respect to nstart
        n3 = MOD((i3+nstart-ip),segbnd) + 1   !next position respect to nstart

        CALL new_pline(pln)   !allocate a new boundary line
        pln%nfirst = nfirst   !associate first position
        nlast  = i1       ! i2-1 = i1       !last point with the same restriction
        IF   (ANY(ifx_b(1:nd1,n1) /= 0) .AND. &
              ANY(ifx_b(1:nd1,n2) /= 0))THEN
          pln%ifix(1:nd1) = 0 !initialize
          DO i = 1,nd1  !assign only if line have same restrictions in nodes
            IF(ifx_b(i,n1) == ifx_b(i,n2)) pln%ifix(i) = ifx_b(i,n1)
          END DO
        ELSE
          pln%ifix(1:nd1) = 0
        END IF
        pln%nlast = nlast                  !last node in the line
        pln%fside = (fside(n1) == 1)    !.TRUE. if side es a free contour
        nplins = nplins + 1                 !increase number of lines
        CALL add_pline(pln,headp,tailp)     !add to the list
        nfirst = nlast + 1                  !new first node for next line
        i1 = nfirst                         !new previous node
        i2 = i1 + 1                         !new present node
        i3 = i2 + 1                         !new next node
        IF (i3 > segbnd) i3 = frstn     !if last node in the list set segbnd
        IF (nlast == segbnd) EXIT       !all nodes processed
      END DO

    END IF

    frstn = frstn + nlbnd(l)
    frsts = frsts + nlbnd(l)
  END DO !end loop over each line in boundary

RETURN
END SUBROUTINE plines
