SUBROUTINE dedge0(nedge,ndime,iwrit,heade,taile)
!     Reads line distributed loads
USE c_input,ONLY: ludat, lures, listen, param, exists, getint, rdfrre, backs, nwopa
USE loa_db,ONLY: edg_nod, new_edg, add_edg
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: nedge, ndime,  iwrit
  TYPE(edg_nod),POINTER:: heade, taile
  !--- Constants
  INTEGER(kind=4),PARAMETER:: nnmax=6
  !--- Local variables
  INTEGER(kind=4):: i, nn, kount, noprs(nnmax), nnode
  REAL(kind=8):: press(ndime,nnmax), xa(ndime)
  TYPE(edg_nod),POINTER:: edg
  CHARACTER(len=1):: systm, defsys
  LOGICAL:: flag=.FALSE.

  !*** distributed edge loads
  CALL listen ('DEDGE0')

  nedge = 0
  IF ( .NOT.exists('EDGE  ')) THEN
    backs = .TRUE.
  ELSE
    IF (iwrit==1) WRITE(lures,"(//5x,'List of Loaded Edges and Applied Loads'/)",ERR=9999)
    IF( exists('NNODE ',i))THEN
      nnode = INT(param(i))
    ELSE
      nnode = 2
    END IF
    IF (exists('L     ')) THEN
      defsys = 'L'
    ELSE IF (exists('G     ')) THEN
      defsys = 'G'
    ELSE
      defsys = 'L'
    END IF

    !*** loop over each loaded edge
    NULLIFY(heade,taile)
    DO
      CALL listen ('DEDGE0')
      IF (exists('ENDEDG')) EXIT

      !*** READ DATA locating the loaded edge and applied load
      nedge = nedge + 1
      IF (exists('L     ',i)) THEN
        systm = 'L'
      ELSE IF (exists('G     ',i)) THEN
        systm = 'G'
      ELSE
        systm = defsys
        i = 0
      END IF

      IF (exists('N1    ',i)) THEN
        noprs(1) = INT(param(i))
        IF (exists('N2    ',i)) THEN
          noprs(2) = INT(param(i))
          nn = 2
        ELSE
          CALL runend('DEDGE0: Specifiy Second node N2')
        END IF
        IF (exists('N3    ',i)) THEN
          noprs(3) = INT(param(i))
          nn = 3
        END IF
      ELSE
        nn = nnode
        noprs(1:nn) = INT(param(i+1:i+nn))
        i = i+nn
      END IF
      xa = 0d0
      IF (ndime == 3 .AND. systm == 'L') THEN
        IF (exists('X     ',i)) xa(1)=param(i)
        IF (exists('Y     ',i)) xa(2)=param(i)
        IF (exists('Z     ',i)) xa(3)=param(i)
      END IF
      IF( nwopa > i )THEN  !more arguments in the line
        kount = nwopa - i
        CALL vecasi(kount,param(i+1),press)
      ELSE                 !read pressure in the next line
        CALL rdfrre('DEDGE0',press,kount,ndime*nn,flag)
      END IF
      IF (kount == ndime) THEN     !generate if only first node is read
        DO i=2,nn
          press(1:ndime,i) = press(1:ndime,1)
        END DO
      END IF
      IF (iwrit == 1) THEN
        IF(nn == 2)WRITE(lures,"(5X,' nodes ',2I7,6F10.3)",ERR=9999)noprs(1:2),press(1:ndime,1:2)
        IF(nn == 3)WRITE(lures,"(5X,' nodes ',3I7,9F10.3)",ERR=9999)noprs(1:3),press(1:ndime,1:3)
      END IF

      CALL new_edg(edg)
      edg%systm2 = systm
      edg%nn2 = nn
      edg%press2(1:ndime,1:nn) = press(1:ndime,1:nn)
      edg%lnod2(1:nn) = noprs(1:nn)
      edg%xa(1:ndime) = xa(1:ndime)
      CALL add_edg(edg, heade, taile)
    END DO
  ENDIF

RETURN
 9999 CALL runen2('')
END SUBROUTINE dedge0
