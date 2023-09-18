SUBROUTINE inidb0a (npoin, nelem, coora, iwork, ndbpoin, nrt,     &
                    icode, xl, head, dtol, lnd)
  !--------------------------------------------------------------
  !gets drawbead line points => LND and assign xy-coords XL 
  !find element to be considered => ICODE  & NRT 
  !--------------------------------------------------------------

  !   m o d u l e s

  USE surf_db

  IMPLICIT NONE

  !   d u m m y   a r g u m e n t s

  INTEGER, INTENT (IN) :: npoin, & !total nodes in the model
                          nelem, & !number of sheet elements
                          ndbpoin  !number of nodes in the db line
  INTEGER, INTENT (OUT) :: nrt     !number of elements to consider
  INTEGER :: iwork (:)             !auxiliar
  INTEGER (kind=4), INTENT (OUT) :: lnd (:,:) !db line connectivities
  REAL (kind=8), INTENT(IN) :: dtol           !maximum distance
  REAL (kind=8), INTENT (IN)  :: coora (:,:)  !present coordinates
  REAL (kind=8), INTENT (OUT)  :: xl (:,:)    !drawbead line coord.
  INTEGER, POINTER :: icode (:)     !(nelem) 1:check element 0:do not
  TYPE (srf_seg), POINTER :: head             !pointer to first segment

  !-----------------------------------------------

  !   l o c a l   v a r i a b l e s

  INTEGER :: i, j, chnode
  REAL (kind=8) :: xp, yp, dmin
  TYPE (srf_seg), POINTER :: seg
  !-----------------------------------------------

  !  gets drawbead line points and assign xy-coords
  CALL db3lin (ndbpoin,lnd, xl, coora)

  iwork = 0    !initializes auxiliar array
  seg => head  !point to first sheet element
  DO i=1,nelem    !loop over each element
    ! change node labels into internal numbers
    ! and mark used nodes
    DO j = 1,4   !for each node
      IF ( seg%nodes(j) /= 0) THEN
         seg%nodes(j) = chnode (seg%nodes(j))
         iwork (seg%nodes(j)) = 1  !mark used node
      END IF
    END DO
    IF ( seg%nodes(4) == 0) seg%nodes(4) = seg%nodes(3)

    seg => seg%next  !point to next element
  END DO

  DO i = 1,npoin  !for each node in the model
    IF (iwork(i) == 1) CYCLE  !not used
    ! check distance
    xp = coora (1, i)  !x coord
    yp = coora (2, i)  !y coord
    dmin = 1d30
    DO j = 1,ndbpoin   !find minimun distance to DB line
      dmin = MIN (dmin, length (xp,yp,xl(1,j),xl(2,j)) )
    END DO
    IF (dmin > dtol) iwork(i) = 0 ! outside bounding box for the drawbead
  END DO

  ALLOCATE ( icode(nelem) )  !auxiliar array
  icode = 1                  !initializes: check all

  seg => head         !point to first element
  DO i=1,nelem          !loop
    ! select sheet segments inside bounding box
    IF (ANY(iwork(seg%nodes(:)) == 0) ) icode(i) = 0 !do not check

    seg => seg%next   !point to next element
  END DO

  nrt = COUNT (icode(:)==1)  !number of elements to be considered

  RETURN

END SUBROUTINE inidb0a

FUNCTION length (xp,yp,xl,yl)
  REAL (kind=8) :: length,xp,yp,xl,yl

  REAL (kind=8) :: dx,dy

  dx = xp-xl
  dy = yp-yl
  length = dx*dx + dy*dy
  length = SQRT (length)
  RETURN
END FUNCTION length
