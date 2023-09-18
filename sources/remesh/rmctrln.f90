SUBROUTINE rmctrln(nelem,heade,nelbnd,cdbnd,lnbnd,nlbnd,fside)
!-----------------------------------------------------------------------------------
!  Create a conectivities and coordinates array for contour line (closed continuous line)
!  Routine works for 2-node segments only (2D problem)
!-----------------------------------------------------------------------------------
USE npo_db,ONLY: coora,label
USE surf_db,ONLY: srf_seg
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: nelem     !number of elements
  INTEGER(kind=4),INTENT(OUT):: nelbnd   !Number of contour segments
  INTEGER(kind=4),POINTER:: lnbnd(:)     !Conectivities of nodes segments
  INTEGER(kind=4),POINTER:: nlbnd(:)     !Number of Continuos lines
  INTEGER(kind=4),POINTER:: fside(:)     !Number of Continuos lines
  REAL(kind=8),POINTER:: cdbnd(:,:)      !Coordinates of nodes segments
  TYPE(srf_seg),POINTER:: heade
  !--- Constants
  INTEGER(kind=4),PARAMETER:: ndime=2
  !--- Local variables
  INTEGER(kind=4):: i, j, k, l, nxnpt, frstn, nline
  TYPE(srf_seg),POINTER:: seg

  !Determination of the total number of continuous lines in the boundary
  !at THIS MOMENT the segments must be in order ------------  WC
  seg => heade       ! point to first segment
  nline = 0             ! number of continuos lines
  frstn = seg%nodes(1)      ! first node in a continuos line
  DO i=1,nelem
    IF( .NOT.ASSOCIATED(seg) )EXIT
    nxnpt = seg%nodes(2)      !label to find the end of line
    IF (nxnpt == frstn ) THEN
      nline = nline + 1  !exist a new continuos line
      IF( ASSOCIATED(seg%next) )THEN
        seg => seg%next        !next segment in contour
        frstn = seg%nodes(1) !new first node in a line
        CYCLE  !find a new continuos line
      END IF
    END IF
    seg => seg%next !next segment in the list
  END DO

  IF( nline == 0 ) CALL runen3('RMCTRLN: CERO LINES IN CONTOURS')

  nelbnd = nelem
  ALLOCATE(cdbnd(ndime,nelbnd),lnbnd(nelbnd),nlbnd(nline),fside(nelbnd))
  cdbnd = 0d0
  nlbnd = 0
  lnbnd = 0
  fside = 0

  !coordinates evaluation
  seg => heade              ! point to first segment
  l = 1                     ! number of line
  k = 0                     ! label nodes position counter
  frstn = seg%nodes(1)      ! first node in a continuos line
  DO i=1,nelem              ! this seems to work for 2-node segments only that is to say 2-D problems only
    cdbnd(1:ndime,i) = coora(1:ndime,seg%nodes(1)) ! used later to analyze contour vertices
    nlbnd(l) = nlbnd(l) + 1     !add a new node to continuos line
    lnbnd(i) = seg%nodes(1)     !Stampack node label
    IF(seg%frees) fside(i) = 1  !is a free side

    nxnpt = seg%nodes(2)       !label to find next segment
    IF (ASSOCIATED(seg%next)) seg=>seg%next      !next segment in the list
    IF (nxnpt == seg%nodes(1)) CYCLE  !IF to ensure continuous line
    IF ((nxnpt == frstn).AND.(l < nline)) THEN ! IF found a line end
      l = l + 1              ! increment line counter and find new segments
      frstn = seg%nodes(1)   ! new node for line start
      CYCLE
    END IF
    !Find next continuous segment from the list
    seg => heade
    DO j=1,nelem
      IF (nxnpt == seg%nodes(1)) EXIT   !If next element is found exit
      seg => seg%next
    END DO
  END DO

RETURN
END SUBROUTINE rmctrln
