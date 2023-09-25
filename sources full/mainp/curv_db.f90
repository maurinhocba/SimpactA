 MODULE curv_db

   ! stores information of curves used for Loads & Prescribed Velocities
   ! created/modified in RDLOAD & RDVEL0
   ! used by APLFOL EXPLIT
   ! by DUMPIN RESTAR

   USE param_db,ONLY: mnam
   IMPLICIT NONE

   SAVE

   TYPE curpar    !Definition: List of nodes (edges, level structure, etc)
     CHARACTER (len=mnam) :: name ! curve name (label)
     INTEGER(kind=4):: dim
     REAL(kind=8),POINTER:: vec(:)
     TYPE (curpar),POINTER:: next
   END TYPE curpar

   INTEGER (kind=4) ::  nlcur   !number of curves
   INTEGER (kind=4) ::  numpac  !Total number of parameters defining curves
   CHARACTER (len=mnam), POINTER :: curnam(:) !(numpac) curve names (labels)
   REAL (kind=8),POINTER :: curpa(:)          !(numpac) parameters defining curves
   TYPE(curpar),POINTER:: headcp, tailcp     !pointers to curve data

 CONTAINS

   FUNCTION getcun (cnam)

     !   returns curve number

     USE lispa0, ONLY :  lures
     IMPLICIT NONE
     INTEGER (kind=4)  :: getcun               ! curve number
     CHARACTER (len=*), INTENT(IN) :: cnam     ! curve name

     LOGICAL          :: found
     INTEGER (kind=4) :: i

     found = .FALSE.
     DO i = 1,nlcur
       found = (TRIM(cnam) == TRIM(curnam(i))) ! compare searched name with the name in db
       IF (found) THEN
         getcun = i                ! returned number of the curve found
         RETURN
       END IF
     END DO

     WRITE(lures,"('Curve label',A)",ERR=9999) TRIM(cnam)
     WRITE(*    ,"('Curve label',A)") TRIM(cnam)
     CALL runend('GETCUN: Curve Label not found')

   RETURN
   9999 CALL runen2('')
   END FUNCTION getcun

   SUBROUTINE dump_curv

     !   dumps curve database
     IMPLICIT NONE

     WRITE (50,ERR=9999) nlcur, numpac   !number of curves  & size of curpa
     WRITE (50,ERR=9999) curnam, curpa   !curve names & curve parameters

     RETURN
     9999 CALL runen2('')
   END SUBROUTINE dump_curv

   SUBROUTINE rest_curv

     !   restores curve database

     IMPLICIT NONE
     TYPE(curpar),POINTER:: aux
     INTEGER :: i,ib,ie

     READ (51) nlcur, numpac     !number of curves  & size of curpa
     ALLOCATE ( curnam(nlcur), curpa(numpac) ) !get memory
     READ (51) curnam, curpa     !curve names & curve parameters
     ! Construct data base from data
     DO i=1,nlcur                !for each curve
       ALLOCATE (aux)            !get memory
       aux%name= curnam(i)       !pass curve name
       ib = curpa(i)             !first position in curpa
       IF( i == nlcur )THEN      !for the last curve
         ie = numpac             !last position in curpa
       ELSE
         ie = curpa(i+1) - 1     !last position for this curve
       END IF
       aux%dim = ie - ib + 1     !number of parameters for this curve
       ALLOCATE( aux%vec(aux%dim) )   !get memory
       aux%vec = curpa(ib:ie)    !pass curve parameters
       CALL add_cur(aux,headcp,tailcp)   !add in list
     END DO

     RETURN

   END SUBROUTINE rest_curv

   SUBROUTINE gen_curpa

     !   generates array curpa

     IMPLICIT NONE
     TYPE(curpar),POINTER:: aux
     INTEGER :: i,ib,ie

     IF( ASSOCIATED(curnam) )DEALLOCATE(curnam) !remove if array exists
     IF( ASSOCIATED(curpa ) )DEALLOCATE(curpa ) !remove if array exists
     nlcur  = 0          !initializes number of curves
     numpac = 0          !initializes number of parameters
     aux => headcp       !point to first curve
     DO              !first loop for each curve
       IF( .NOT. ASSOCIATED(aux) )EXIT  !exit if last curve
       nlcur  = nlcur  + 1              !increase number of curves
       numpac = numpac + aux%dim + 1    !increase number of parameters
       aux => aux%next                  !point to next curve
     END DO

     ALLOCATE (curnam(nlcur), curpa(numpac) ) !get memory
     aux => headcp       !point to first curve
     ib = nlcur+1        !pointer to first curva parameters
     DO i=1,nlcur    !second loop for each curve
       curnam(i) = aux%name     !curname
       curpa(i)  = ib           !pointer to curve parameters in curpa
       ie = ib + aux%dim - 1    !last position
       curpa(ib:ie) = aux%vec   !pass arguments
       IF( INT(curpa(ib)) == 3 )curpa(ib+4) = ib+5 !initializes actual position
       ib = ie + 1              !new first position
       aux => aux%next          !point to next curve
     END DO

     RETURN

   END SUBROUTINE gen_curpa

!-----------------------------------------------------------------------
  SUBROUTINE add_cur(new,head,tail)
!-----------------------------------------------------------------------
  !This subroutine adds data to the end of the list
  TYPE(curpar),POINTER:: new, head, tail

  IF (.NOT.ASSOCIATED(head)) THEN
    !list is empty, start it
    head => new
    tail => new
  ELSE
    !Check if a list is empty
    tail%next => new
    tail => new
  END IF
  NULLIFY (tail%next)

  RETURN
  END SUBROUTINE add_cur

!-----------------------------------------------------------------------
  SUBROUTINE del_cur(curname)
!-----------------------------------------------------------------------
  !Delete a curve from the list
  IMPLICIT NONE

  CHARACTER (len=mnam), INTENT(IN) :: curname     ! curve name
  TYPE(curpar),POINTER :: aux,prev

  aux => headcp
  NULLIFY(prev) !first curve in the list
  DO
    IF( .NOT.ASSOCIATED (aux) )THEN
      WRITE(55,"(' curve ',A,' does not exists ')",ERR=9999) TRIM(curname)
      RETURN
    END IF
    IF( TRIM(curname) == TRIM(aux%name) )THEN
      IF( ASSOCIATED(prev)) THEN  !not the first curve in the list
        prev%next => aux%next    !
        IF( .NOT.ASSOCIATED (prev%next) )tailcp => prev
      ELSE
        headcp => aux%next
      END IF
      DEALLOCATE(aux%vec)
      DEALLOCATE(aux)
      EXIT
    END IF
    prev => aux
    aux => aux%next
  END DO

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE del_cur

!!-----------------------------------------------------------------------
!  SUBROUTINE dalloc_cur(head,tail)
!!-----------------------------------------------------------------------
!  !Delete a list
!  TYPE(curpar),POINTER:: head,tail
!  !Local variables
!  TYPE(curpar),POINTER:: aux
!
!  IF (.NOT.ASSOCIATED(head)) RETURN
!
!  !Search the starting point
!  NULLIFY(tail)
!  DO
!    IF (.NOT.ASSOCIATED(head)) EXIT
!    aux => head%next
!    CALL del_cur(head)
!    head => aux
!  END DO
!
!  RETURN
!  END SUBROUTINE dalloc_cur

 END MODULE curv_db
