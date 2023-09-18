SUBROUTINE rnmesh(scalef,f_axsym)
!=> remesh_db   coorn(ndime,numnp)  l_new(nnode,ne_new)
! READ new mesh: nodal coordinates and element connectivities
USE ctrl_db, ONLY: ndime
USE c_input
USE meshmo_db,ONLY: numpn, ne_new, l_new, coorn, r_elm_zone
IMPLICIT NONE

  !--- Dummy variables
  REAL(kind=8),INTENT(IN):: scalef
  LOGICAL,OPTIONAL:: f_axsym  ! TRUE => axisym. problem, introduced to correct negative x
  !--- Local variables
  INTEGER(kind=4):: n, nn, ln(4), iwrit=0,  nnode
  REAL(kind=8):: area, det2, xe(ndime,3)
  ! new variables
  INTEGER(kind=4), ALLOCATABLE :: lc(:,:)
  INTEGER(kind=4) :: i,j,k,l,ei,ej,ii,ij,jj
  LOGICAL :: found
  INTEGER (kind=4), PARAMETER :: nextn(3) = (/2,3,1/) !cycling list

  CALL listen('RNMESH')                !read first control card
  numpn = getint('NUMNP ',  0,'!Number of new points  ........... ')
  ALLOCATE(coorn(ndime,numpn))         !new coordinates array
  IF(iwrit == 1) WRITE(lures,"(//'   ECHO of the Coordinate Data ',/ &
  & 5x,' Node      X',9X,'Y',9X,'Z')",ERR=9999)

  ! nodes are read in sequential order, from 1 to NUMPN
  DO n=1,numpn                      !for each new node
    CALL listen('RNMESH')                          !read a card
    coorn(1:ndime,n) = param(2:ndime+1)/scalef     !transfer coordinates
    IF (f_axsym) coorn(1,n) = MAX(coorn(1,n),0d0)  ! correct neg.coord
    IF (iwrit==1) WRITE(lures,"(i10,3f10.5)",ERR=9999) n,param(2:ndime+1)  !ECHO
  END DO

  CALL listen('RNMESH')                !read second control card
  nnode = getint('NNODE ',  0,'!Number of nodes per element ..... ')
  ne_new = getint('NELEM ',  0,'!Number of elements in a new mesh  ')
  ALLOCATE(l_new(nnode,ne_new))        !new connectivites
  ! elements are read in sequential order, from 1 to ne_new
  DO n=1,ne_new                !for each new element
    CALL listen('RNMESH')                    !read a card
    !    note that associated material must not be present
    IF (nnode == 4) THEN
      l_new(1:nnode,n) = param(2:nnode+1)     !transfer connectivities
    ELSE ! IF (nnode == 3)
      ln(1:nnode) = param(2:nnode+1)          !transfer connectivities
      ! check normals - GiD 6 gives unstable results
      ! after GiD is corrected this part can be removed
      xe(:,:) = coorn(:,ln(1:nnode))          !triangle coordinates
      area  = det2(xe(1,1),xe(1,2),xe(1,3),ndime,.FALSE.)
      IF (area < 0) THEN
        ! clockwise numeration swapped
        DO nn=1,nnode
          l_new(nnode+1-nn,n) = ln(nn)     !transfer connectivities
        END DO
      ELSE
        l_new(1:nnode,n) = ln(1:nnode)
      END IF
    END IF
  END DO

  !  swap corner elements for triangles (not in zone remeshing)
  IF( nnode == 3 .AND. .NOT.r_elm_zone ) THEN
    ALLOCATE(lc(3,ne_new))        !new connectivites
    lc = 0
    DO ei=1,ne_new
      DO ii=1,3                        !for each side in the element
        IF( lc(ii,ei) /= 0 )  CYCLE    !if opposite already find, CYCLE
        ij = nextn(ii)                 !next local node
        j  = l_new(ij,ei)              !first node of the side
        ij = nextn(ij)                 !next local node
        k  = l_new(ij,ei)              !second node of the side
        found = .FALSE.                !initializes search
        DO ej=ei+1,ne_new   !search an element with the same side (node1 -- node2)
          DO jj=1,3                         !for each side
            IF( lc(jj,ej) /= 0 )CYCLE  !if side already paired, CYCLE
            ij = nextn(jj)             !next local node
            IF(l_new(ij,ej) /= k)CYCLE    !first test, if failed, CYCLE
            ij = nextn(ij)             !next local node
            IF(j == l_new(ij,ej)) THEN            !second test
              found = .TRUE.           !side found
              lc(ii,ei) = l_new(jj,ej) !opposite global node
              lc(jj,ej) = l_new(ii,ei) !opposite global node
            END IF
            EXIT                            !found or not exit loop
          END DO
          ! if twin side found  EXIT loop
          IF (found) EXIT
        END DO
      END DO
    END DO
    ! check elements with only one neighbour
    DO ei=1,ne_new
      nn = 0
      DO ii=1,3                        !for each side in the element
        IF( lc(ii,ei) /= 0 ) nn = nn+1
      END DO
      IF( nn == 1) THEN         !element has one neighbour only
        ii = 1                           !initializes loop counter
        DO                                  !for each side in the element
          IF( lc(ii,ei) /= 0 )EXIT   !side found, EXIT
          ii = ii + 1                 !increase counter
        END DO
        i  = l_new(ii,ei)                   !only node with opposite element
        ij = nextn(ii)                !next local node in element I
        j  = l_new(ij,ei)             !first node of the side I
        ij = nextn(ij)                !next local node
        k  = l_new(ij,ei)             !second node of the side I

        loop1 : DO ej=1,ne_new   !search an element with node i as the opposite node
          jj = 1                            !initializes loop counter
          DO                                !for each side in the element
            IF( jj > 3 )EXIT                !number of nodes exceeded, EXIT inner loop
            IF( lc(jj,ej) == i  )EXIT loop1 !node found, EXIT outer loop
            jj = jj + 1                    !increase side counter
          END DO
        END DO loop1
        IF( ANY(lc(:,ej) == 0))CYCLE  !check if the other element has 3 neighbours
        !swap elements
        l  = lc(ii,ei)                !extra node (opposite to i)
        l_new(:,ei) = (/ i,j,l /)
        l_new(:,ej) = (/ k,i,l /)
      END IF
    END DO

  END IF
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE rnmesh
