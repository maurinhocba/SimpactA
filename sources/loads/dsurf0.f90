SUBROUTINE dsurf0(nsurf,ndime,iwrit,heads,tails)
!Reads surface distributed loads
USE c_input,ONLY: ludat, lures, listen, param, words, exists, rdfrre, backs, nwopa
USE loa_db,ONLY: srf_nod, new_srf, add_srf
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: nsurf, ndime, iwrit
  TYPE(srf_nod),POINTER:: heads, tails
  !--- Constants
  INTEGER(kind=4),PARAMETER:: nnmax=6
  !--- Local variables
  INTEGER(kind=4):: i, n, nn, nvert, ivert, kount, noprs(nnmax), nnode
  REAL(kind=8):: direc(ndime), p(nnmax), area
  TYPE(srf_nod),POINTER:: srf
  CHARACTER(len=1):: systm, defsys
  LOGICAL:: flag=.FALSE.

  !*** distributed surface loads
  CALL listen('DSURF0')

  nsurf = 0

  IF( .NOT.exists('SURF  ')) THEN
    backs = .TRUE.
  ELSE

    IF(iwrit == 1) WRITE(lures,"(/6x,'List of loaded Faces read in this strategy'/)",ERR=9999)
    IF( exists('NNODE ',i)) THEN
      nnode = INT(param(i))
    ELSE
      nnode = 3
    END IF
    IF (exists('L     ')) THEN
      defsys = 'L'
    ELSE IF (exists('G     ')) THEN
      defsys = 'G'
    ELSE
      defsys = 'L'
    END IF
    !*** loop over each loaded face
    NULLIFY(heads,tails)
    p = 0d0    !initializes
    noprs = 0
    DO
      CALL listen('DSURF0')
      IF (exists('ENDSUR')) EXIT

      nsurf = nsurf + 1
      IF (exists('L     ',i)) THEN
        systm = 'L'
      ELSE IF (exists('G     ',i)) THEN
        systm = 'G'
      ELSE
        systm = defsys
        i = 0
      END IF
      IF (exists('NN    ',i)) THEN
        nn = INT(param(i))
      ELSE
        nn = nnode
      END IF
      DO n=1,nn
        IF (LEN_TRIM(words(i+n)) == 0) THEN
          noprs(n) = INT(param(i+n))
        ELSE
          CALL runend('DSURF0: ERROR reading face nodes')
        END IF
      END DO
      i = i + nn
      nvert = nn
      IF( nn > 4) nvert=nn/2

      direc = 0d0
      IF (exists('DX    ',i)) direc(1)=param(i)
      IF (exists('DY    ',i)) direc(2)=param(i)
      IF (exists('DZ    ',i)) direc(3)=param(i)
      CALL vecuni(ndime,direc,area)
      IF (area == 0d0) direc(3)=1d0
      IF( nwopa > i )THEN  !more arguments in the line
        kount = nwopa - i
        p(1:kount) = param(i+1:nwopa)
      ELSE                 !read pressure in the next line
        CALL rdfrre('DSURF0',p,kount,nvert,flag)
      END IF
      IF (kount == 1) THEN
        p(2:nn) = p(1)
      ELSE IF (nn > nvert) THEN
        DO ivert=nvert+1,nn-1
          p(ivert) = (p(ivert-nvert)+p(ivert-nvert+1))/2d0
        END DO
        p(nn) = (p(1)+p(nvert))/2d0
      END IF

      !prints actual DATA
      IF (iwrit == 1) THEN
        WRITE(lures,"(' FACE',i3,' Sys= ',a1,' NN=',i2,' DIR',3f6.3,' Nodes',6i6)",ERR=9999)   &
           nsurf,systm,nn,direc,noprs(1:NN)
        WRITE(lures,"(' Nodal p ',6g12.3)",ERR=9999) p(1:nn)
      END IF

      CALL new_srf(srf)
      srf%systm3 = systm
      srf%nn3 = nn
      srf%press3 = p
      srf%lnod3 = noprs
      srf%direc3(1:ndime) = direc(1:ndime)
      CALL add_srf(srf, heads, tails)
    END DO
  END IF

RETURN
9999 CALL runen2('')
END SUBROUTINE dsurf0
