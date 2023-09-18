MODULE Form_Patch_Lib

  ! library of routines to form a patch of elements

  USE SPR_Elem_Type_Db

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Count_L_at_Node (nelem, nnode, npoin, nodes, l_to_n_sum)

    ! count number of elements adjacent to each node

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: &
      nelem,             & ! number of elements
      nnode,             & ! number of nodes per element
      npoin,             & ! number of nodal points
      nodes (nelem, nnode) ! system array of incidences of all elements
    INTEGER, INTENT(OUT) :: &
      l_to_n_sum (npoin)   ! number of elem neighbors of node

    ! local
    INTEGER :: elem_nodes (nnode)
    INTEGER :: ie, in, n_temp

    l_to_n_sum = 0             ! initialize
    DO ie = 1, nelem           ! loop over all elements
      elem_nodes = get_elem_nodes (ie, nnode, nodes)
      DO in = 1, nnode         ! loop over each node
        n_temp = elem_nodes (in)
        l_to_n_sum (n_temp) = l_to_n_sum (n_temp) + 1
      END DO ! over in
    END DO ! over elements

  END SUBROUTINE Count_L_at_Node


  SUBROUTINE Count_Elems_at_Elem (nelem, nnode, npoin, l_first, l_last, &
                                  nodes, needs, l_to_l_sum)

    ! count number of elements adjacent to other elements

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: &
      nelem,             & ! total number of elements
      nnode,             & ! number of nodes per element
      npoin,             & ! total number of nodes
      needs,             & ! number of common nodes to be a neighbor
      l_first (npoin),   & ! element where node first appears
      l_last  (npoin),   & ! element where node last appears
      nodes (nelem, nnode) ! system array of incidences of all elements
    INTEGER, INTENT(OUT) :: &
      l_to_l_sum (nelem)   ! number of elem neighbors

    ! local
    INTEGER :: found, ie, in, l_test, l_start, l_stop, n_test, &
               nulls, need, kount
    INTEGER :: elem_nodes (nnode), neig_nodes (nnode)

    need       = MAX (1, needs)   ! initialize
    l_to_l_sum = 0

    Main : DO ie = 1, nelem      ! loop over all elements
      found = 0

      ! extract incidences list for element ie
      elem_nodes = get_elem_nodes (ie, nnode, nodes)

      ! establish range of possible element neighbors
      l_start = nelem ; l_stop = 0
      DO in = 1, nnode
        n_test  = elem_nodes (in) !b
        IF ( n_test < 1 ) CYCLE! to a real node
        l_start = MIN (l_start, l_first (n_test) )          !b
        l_stop  = MAX (l_stop,  l_last  (n_test) )          !b
      END DO

      ! loop over possible element neighbors
      IF ( l_start <= l_stop) THEN
        Range : DO l_test = l_start, l_stop
          IF ( l_test /= ie) THEN
            kount = 0  ! no common nodes

            ! loop over incidences of possible element neighbor
            neig_nodes = get_elem_nodes (l_test, nnode, nodes)
            Local : DO in = 1, nnode
              n_test = neig_nodes (in)
              IF ( l_first (n_test)  > ie ) CYCLE Local ! to next node
              IF ( l_last  (n_test)  < ie ) CYCLE Local ! to next node

              ! compare with incidences of element ie
              IF ( ANY ( elem_nodes == n_test ) ) THEN
                kount = kount + 1 ! number of common nodes
                IF ( kount == need ) THEN ! is a neighbor
                  found = found + 1
                  EXIT Local ! this l_test element search loop
                END IF ! number needed
              END IF
            END DO Local ! over in
          END IF
        END DO Range ! over candidate element l_test
      END IF ! a possible candidate
      l_to_l_sum (ie) = found
    END DO main ! over all elements

    !PRINT *, 'Note: maximum number of element neighbors = ', &
    !           MAXVAL ( l_to_l_sum )
    nulls = COUNT ( l_to_l_sum == 0 ) !  check data
    IF ( nulls > 0 ) THEN
      PRINT *, 'Warning, ', nulls, ' elements have no element neighbors'
    END IF
  END SUBROUTINE count_elems_at_elem

  SUBROUTINE First_Last_L_at_Pt (nodes, l_first, l_last)

  ! find element of first & last appearance of each node

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: &
    nodes (nnode, nelem)    ! systen array of all element incidences
  INTEGER, INTENT(OUT) :: &
    l_first (npoin),      & ! element of first appearance of node i
    l_last (npoin)          ! element of last appearance of node i

  INTEGER, ALLOCATABLE :: el_nodes (:)
  INTEGER :: i, ie, l

   l_first = 0 ; l_last = 0  ! initialize
   ALLOCATE ( el_nodes (nnode) )

   DO ie = 1, nelem

     el_nodes = get_elem_nodes (ie, nnode, nodes)

     DO i = 1, nnode
       l = el_nodes (i)
       IF ( l < 1 ) CYCLE ! to next node
       IF ( l_first (l) == 0 ) l_first (l) = ie
       l_last (l) = ie
     END DO
   END DO

   DEALLOCATE ( el_nodes )

  END SUBROUTINE First_Last_L_at_Pt

  SUBROUTINE Form_Elems_at_El (nelem, nnode, npoin, l_first, l_last,    &
                               nodes, l_to_l_sum, l_to_l_neigh,  &
                               neigh_l, needs)

    ! form list of elements neighbouring to other elements

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: &
      nelem,   & ! total number of elements
      nnode,   & ! number of nodes per element
      npoin,   & ! total number of nodes
      neigh_l, & ! maximum number of neighbors at an element
      needs,   & ! number of common nodes to be a neighbor
      l_first    (npoin),        & ! element where node first appears
      l_last     (npoin),        & ! element where node last appears
      nodes      (nelem, nnode), & ! system array of incidences of all elements
      l_to_l_sum (nelem)           ! number of elem neighbors of element
    INTEGER, INTENT(OUT)   :: &
      l_to_l_neigh (neigh_l, nelem)! elem neighbor j of element i

    ! local
    INTEGER :: ie, in, l_test, l_start, l_stop, n_test, sum_l_to_l, &
      found, next, & ! add need for face or edge sub-sets
      kount, need, where, & !b
      elem_nodes (nnode), neig_nodes (nnode)

    need         = MAX (1, needs)   ! initialize
    l_to_l_neigh = 0                ! initialize
    Main : DO ie = 1, nelem         ! loop over all elements

      sum_l_to_l = l_to_l_sum (ie)               ! max possible neighbors
      found = COUNT ( l_to_l_neigh (:, ie) > 0 ) ! previously found
      IF ( found == sum_l_to_l ) CYCLE main      ! all have been found

      elem_nodes = get_elem_nodes (ie, nnode, nodes)

      ! establish range of possible element neighbors
      l_start = nelem + 1 ; l_stop = 0
      DO in = 1, nnode
        l_start = MIN (l_start, l_first (elem_nodes (in)) )
        l_stop  = MAX (l_stop,  l_last  (elem_nodes (in)) )
      END DO
      l_start = MAX (l_start, ie+1) ! search above ie only

      ! loop over possible element neighbors
      IF ( l_start <= l_stop) THEN
        Range : DO l_test = l_start, l_stop
          kount = 0  ! no common nodes

          ! extract nodes of l_test
          neig_nodes = get_elem_nodes (l_test, nnode, nodes)

          ! loop over incidences of possible element neighbor
          Local : DO in = 1, nnode
            n_test = neig_nodes (in)
            IF ( n_test <  1 ) CYCLE ! to a real node
            IF ( l_first (n_test)  > ie ) CYCLE Local ! to next node
            IF ( l_last  (n_test)  < ie ) CYCLE Local ! to next node

            ! compare with incidences of element ie
            IF ( ANY ( elem_nodes == n_test ) ) THEN
              kount = kount + 1         ! shared node count
              IF ( kount == need ) THEN ! neighbor pair found
                found = found + 1       ! insert the pair

                ! note: this insert is not ordered. if a specific face
                !       order is required determine where here

                where = found ! or order the current face
                l_to_l_neigh (where, ie) = l_test     ! one of two

                next = COUNT ( l_to_l_neigh(:, l_test) > 0 ) !b /= ??
                where = next+1 ! or order the neighbor face

                IF ( l_to_l_sum (l_test) > next )  &
                  l_to_l_neigh (next+1, l_test) = ie  ! two of two
                IF ( sum_l_to_l == found ) CYCLE main ! all found
                CYCLE Range ! this l_test element search loop
              END IF ! number needed
            END IF ! share at least one common node
          END DO Local ! over n_test
        END DO Range ! over candidate element l_test
      END IF ! a possible candidate
    END DO main! over all elements

  END SUBROUTINE Form_Elems_at_El

  SUBROUTINE Calc_GP_Coord (nodes, x, coogp, nelem, nnode, npoin, ndime, ngaus)

    ! calculate cartesian coordinates of Gauss points

    IMPLICIT NONE
    INTEGER, INTENT (IN) :: &
      nelem, & ! total number of elements
      nnode, & ! number of nodes per element
      npoin, & ! total number of nodes
      ndime, & ! dimensions size of the problem
      ngaus, & ! number of Gauss points per element
      nodes (nnode, nelem)  ! nodal incidences of all elements
    REAL(kind=8), INTENT (IN ) :: &
      x     (ndime, npoin)        ! coordinates of system nodes
    REAL(kind=8), INTENT (OUT) :: &
      coogp (ndime, ngaus, nelem) ! coordinates of Gauss points

    ! local
    INTEGER  :: ie, ip
    REAL(kind=8) :: xyz (ndime)

    DO ie = 1, nelem ! loop over elements

      ! get element node numbers & coordinates
      !elem_nodes = get_elem_nodes (ie, nnode, nodes)
      elem_nodes (1:nnode) = nodes (1:nnode, ie)
      coord (1:nnode, 1:ndime) = TRANSPOSE( x (1:ndime, elem_nodes (1:nnode)) )
      !call elem_coord (nnode, ndime, x, coord, elem_nodes)

      DO ip = 1, ngaus           ! loop over gauss points
        h   = get_h_at_qp (ip)   ! evaluate interpolation functions
        xyz = matmul (h, coord)  ! find global coord, (isoparametric)
        coogp(:,ip,ie)=xyz       ! store coords
      END DO

    END DO                       ! end of loop over elements

    RETURN

  END SUBROUTINE Calc_GP_Coord

  SUBROUTINE Find_GP_Lproj(l_old, coora, cgp_new, l_proj, ne_new, nnode, ngaus )

    ! finds projection of new Gauss points into old elements

    USE lispa0
    USE ctrl_db, ONLY: ndime, ndofn
    IMPLICIT NONE

    ! dummy arguments
    INTEGER (kind=4) :: &
      ne_new,    &    ! number of new elements
      nnode,     &    ! number of node/element
      ngaus,     &    ! number of node/element
      l_old(:,:)      ! connectivities of old elements
    INTEGER (kind=4) :: &
      l_proj(:,:)      ! element where gp project
    REAL (kind=8)    :: &
      cgp_new(:,:,:), & ! coordinates of new points
      coora(:,:)        ! coordinates of old points

    ! local
    !LOGICAL nonear     ! point not near
    INTEGER (kind=4) :: i, ie, ig, il, ne_old, nepe,  &
                        l(3,2),io
    REAL (kind=8) :: xn(ndime), x(ndime,3), ar(3), area, det2
    REAL (kind=8), PARAMETER ::  tiny=5.0d-2    !tolerance in proyection
    INTERFACE
      INCLUDE 'nonear.h'
    END INTERFACE

    ne_old = SIZE(l_old,2)       !number of elements in previous mesh
    l_proj = 0

    ! loop over new gauss points

    ie = 1                      !initializes first element to look in
    L10: DO il = 1, ne_new      !for each new point
    L11: DO ig = 1, ngaus       !for each new point
      xn  = cgp_new(:,ig,il)    !new coordinates

      ! loop over all the elements in the old mesh to find its proyection
      io = ie                   !remember first searched element

      L15: DO                   !loop until proyection found

        CALL trian(nnode,l_old(1:nnode,ie),l,nepe) !divide element in triangles

        L20: DO i=1,nepe  ! number of sub-triangles in an element

          x(:,:) = coora(:,l(:,i))    !triangle coordinates

          ! discard an element if this point is outside its boundingbox
          IF (nonear (xn, x, ndime)) CYCLE L20

          ! element area
          area  = det2 (x(1,1),x(1,2),x(1,3),ndime,.true.)

          ! three subelements (natural triangular coordinates)
          ar(1) = det2 (x(1,2), x(1,3), xn, ndime,.false.) / area
          IF (ar(1) < -tiny) CYCLE L20

          ar(2) = det2 (x(1,3), x(1,1), xn, ndime,.false.) / area
          IF (ar(2) < -tiny) CYCLE L20

          ar(3) = det2 (x(1,1), x(1,2), xn, ndime,.false.) / area
          IF (ar(3) < -tiny) CYCLE L20

           l_proj(ig,il) = ie
          CYCLE L11  !jump out - don't check any more elements

        END DO L20   !loop over individual triangles
        ie = ie+1                !new element
        IF( ie == io )EXIT       !if back to the first element tested =>exit
        IF( ie > ne_old ) ie = 1 !if passed end of the list, back to first
      END DO L15
      WRITE(lures,*,ERR=9999) ig,' of',il, ' not placed, fatal error'!projection not found
      STOP

    END DO L11  ! end of new point loop
    END DO L10  ! end of new point loop

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE Find_GP_Lproj

END MODULE Form_Patch_Lib

