 MODULE flc_db
 IMPLICIT NONE

   TYPE flc_tp
     INTEGER(kind=4):: lbl,   & !FLC label
                       sttyp, & !strain type (to read & write in postprocess file, see "wrflc.f90")
                                !     1: Unitary strain is written     flc%cv(1:2,i)
                                !     2: Porcentual strain is written  1d2*flc%cv(1:2,i)
                                !     3: True strain is written        LOG(1d0+flc%cv(1:2,i))
                       npt      !Number of points that define FLC curve
     REAL(kind=8),POINTER:: cv(:,:)   !FLC curve defined by points
     REAL(kind=8):: LmM,      & !Marginal tolerance
                    LmPS,     & !Plane strain tolerance
                    LmLS,     & !Low strain tolerance
                    LmT         !Excessive thinning
     TYPE(flc_tp),POINTER:: next
   END TYPE flc_tp

   !--- General variables
   INTEGER(kind=4):: nflc=0
   TYPE(flc_tp),POINTER:: hpflc=>NULL(), tpflc=>NULL()

 CONTAINS

   !------------------------------------------------
   SUBROUTINE ini_flc_db
      !TRUE: Define default FLC (unitary unit)
      TYPE (flc_tp), POINTER :: flc
      REAL(kind=8):: pi

      pi = 4d0*DATAN(1d0)                     !pi
      CALL new_flc(flc)
      flc%npt = 7
      ALLOCATE(flc%cv(2,flc%npt))
      flc%cv(1,1) =-20d-2  ;  flc%cv(2,1) = 53d-2
      flc%cv(1,2) =-10d-2  ;  flc%cv(2,2) = 44d-2
      flc%cv(1,3) =  0d-2  ;  flc%cv(2,3) = 34d-2
      flc%cv(1,4) =  5d-2  ;  flc%cv(2,4) = 28d-2
      flc%cv(1,5) = 10d-2  ;  flc%cv(2,5) = 33d-2
      flc%cv(1,6) = 20d-2  ;  flc%cv(2,6) = 38d-2
      flc%cv(1,7) = 40d-2  ;  flc%cv(2,7) = 47d-2
      flc%LmM = 10d-2
      flc%LmPS = pi/36d0
      flc%LmLS = 5d-3
      flc%LmT = 20d-2
      CALL add_flc(flc,hpflc,tpflc)
      nflc = nflc+1
   END SUBROUTINE ini_flc_db

   !------------------------------------------------
   SUBROUTINE new_flc(new)
   !Create a new 'flc_tp' pointer
   IMPLICIT NONE
   TYPE(flc_tp),POINTER:: new

     ALLOCATE(new)
     new%lbl = 0
     new%sttyp = 1        !     1: Unitary strain is read by default
     new%npt = 0
     new%LmM  = 0d0
     new%LmPS = 0d0
     new%LmLS = 0d0
     new%LmT  = 0d0
     NULLIFY(new%cv,new%next)

   RETURN
   END SUBROUTINE new_flc

   !------------------------------------------------
   SUBROUTINE add_flc(new,head,tail)
   !This subroutine adds data to the end of the list
   IMPLICIT NONE
   TYPE(flc_tp),POINTER:: new, head, tail

     !Check if a list is empty
     IF (.NOT.ASSOCIATED(head)) THEN
       !list is empty, start it
       head => new
       tail => new
     ELSE
       !add a segment to the list
       tail%next => new
       tail => new
     ENDIF

   RETURN
   END SUBROUTINE add_flc

   !------------------------------------------------
   SUBROUTINE del_flc(elm,head,tail)
   !Delete a list
   TYPE(flc_tp),POINTER:: elm, head, tail
   !Local variables
   TYPE(flc_tp),POINTER:: ant, aux

     NULLIFY(ant)
     aux => head
     DO
       IF (.NOT.ASSOCIATED(aux)) EXIT
       IF (aux%lbl == elm%lbl) EXIT
       ant => aux
       aux => aux%next
     END DO

     IF (ASSOCIATED(aux)) THEN
       IF (.NOT.ASSOCIATED(ant)) THEN
         head => elm%next
       ELSE IF (.NOT.ASSOCIATED(elm%next)) THEN
         tail => ant
         NULLIFY(ant%next)
       ELSE
         ant%next => elm%next
       END IF
       NULLIFY(aux)
     END IF

     IF (ASSOCIATED(elm%cv)) DEALLOCATE(elm%cv)
     DEALLOCATE(elm)

   RETURN
   END SUBROUTINE del_flc

 !!  !------------------------------------------------
 !!  SUBROUTINE dalloc_flc(head,tail)
 !!  !Delete a list
 !!  TYPE(flc_tp),POINTER:: head, tail
 !!  !Local variables
 !!  TYPE(flc_tp),POINTER:: aux
 !!
 !!    NULLIFY(tail)
 !!    DO
 !!      IF (.NOT.ASSOCIATED(head)) EXIT
 !!      aux => head%next
 !!      IF (ASSOCIATED(head%cv)) DEALLOCATE(head%cv)
 !!      DEALLOCATE(head)
 !!      head => aux
 !!    END DO
 !!
 !!  RETURN
 !!  END SUBROUTINE dalloc_flc

   !------------------------------------------------
   SUBROUTINE srch_flc(head,lbl,found,flc)
   !This subroutine searches for a FLC curve labed "lbl"
   LOGICAL,INTENT(OUT):: found
   INTEGER(kind=4),INTENT(IN):: lbl
   TYPE(flc_tp),POINTER:: head, flc

     found = .FALSE.
     !Check if a list is empty
     IF (ASSOCIATED (head)) THEN
       flc => head
       DO
         IF (.NOT.ASSOCIATED(flc)) EXIT
         IF (flc%lbl == lbl) THEN
           found = .TRUE.
           EXIT
         END IF
         flc => flc%next
       END DO
     ENDIF
     IF (.NOT.found) NULLIFY(flc)

   RETURN
   END SUBROUTINE srch_flc

   !------------------------------------------------
   SUBROUTINE dump_flc()
   !This subroutine writes data from restart file
   IMPLICIT NONE

   !--- Local variables
   TYPE(flc_tp),POINTER:: flc

   WRITE(50,ERR=9999) nflc

   flc => hpflc
   DO
     IF (.NOT.ASSOCIATED(flc)) EXIT
     WRITE(50,ERR=9999) flc%lbl, flc%sttyp, flc%npt
     IF (flc%npt > 0) WRITE(50,ERR=9999) flc%cv(1:2,1:flc%npt)
     WRITE(50,ERR=9999) flc%LmM, flc%LmPS, flc%LmLS, flc%LmT
     flc => flc%next
   END DO

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE dump_flc

   !------------------------------------------------
   SUBROUTINE rest_flc()
   !This subroutine reads data from restart file
   IMPLICIT NONE

   !--- Local variables
   INTEGER(kind=4):: lbl, sttyp, npt, i
   TYPE(flc_tp),POINTER:: flc

   READ(51) nflc
   DO i=1,nflc
     READ(51) lbl,sttyp,npt
     CALL new_flc(flc)
     flc%lbl = lbl
     flc%sttyp = sttyp
     flc%npt = npt
     IF (flc%npt > 0) THEN
       ALLOCATE(flc%cv(2,flc%npt))
       READ(51) flc%cv(1:2,1:flc%npt)
     END IF
     READ(51) flc%LmM, flc%LmPS, flc%LmLS, flc%LmT
     CALL add_flc(flc,hpflc,tpflc)
   END DO

   RETURN
   END SUBROUTINE rest_flc

 END MODULE flc_db
