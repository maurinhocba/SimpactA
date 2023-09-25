 SUBROUTINE inpd11 (task,iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 2-node 2-d beam/shell element
 !
 !******************************************************************
 USE ele11_db
 USE ctrl_db, ONLY: ntype
 IMPLICIT NONE

  ! dummy arguments
  CHARACTER(len=*):: elsnam    ! element set name
  CHARACTER(len=*):: task      ! requested task
  INTEGER(kind=4)::  nelms,  & ! number of element sets of this type
                     iwrit      ! flag to echo data input

  ! local variables
  LOGICAL:: oldset
  INTEGER(kind=4):: nelem, nstre, nbn, nreqs, narch, ngaus
  CHARACTER(len=mnam):: sname     ! element set name

  CHARACTER(len=12):: ptype(3) =(/ 'Plane stress', &
                                   'Plane strain', &
                                   'Axisymmetric'  /)

  TYPE(ele11_set),POINTER:: elset => null()

  sname = TRIM(elsnam)
  IF (TRIM(task) == 'INPUT') THEN
    ! check if list of sets and set exists and initializes
    CALL srch_ele11 (head, elset, elsnam, oldset)
    IF (oldset) THEN    !if set exists
      CALL comm11 (1, nelem, nstre, nbn, ngaus, nreqs, narch, elsnam, elset)
      elset%lside = .FALSE.  !initializes flag to compute LSIDE
    ELSE                !set ELSNAM does not exist
      nelem = 0
      CALL new_ele11(elset)  !reserve memory for set
      CALL listen('INPD11')  !read a line
      WRITE(lures,"(/,5x,'Control parameters for beam element'// )",ERR=9999)
      ngaus=getint('NGAUS ',2,' Number of Integration points......')
      IF( ngaus == 1 )THEN
        elset%shap = 0.5d0
        elset%stabs = getrea('STABS ',0.15D0,' Stabilization coefficient ........')
      ELSE IF( ngaus == 2 )THEN
        elset%shap(:,1) = (/ 0.788675134d0, 0.211324866d0 /)
        elset%shap(:,2) = (/ 0.211324866d0, 0.788675134d0 /)
        elset%stabs = 0d0
      ELSE
        CALL runend('INPD11, NGAUSS must be 1 or 2 ')
      END IF
      SELECT CASE (ntype)
      CASE (1)
        nstre = 2
      CASE (2:3)
        nstre = 4
      END SELECT
      WRITE(lures,"(10x,A,/,10x,'No of stresses  (NSTRE) =',i3,/)",ERR=9999) &
            TRIM(ptype(ntype)),nstre
      nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
      elset%strai = exists('ZEROCU')
    END IF
    !  read new data or add to previous data
    CALL elmd11(nelem,nstre,ntype,ngaus,elset%head,elset%tail,iwrit)
    IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )
    CALL comm11 (0, nelem, nstre, nbn, ngaus, nreqs, narch, elsnam, elset)
    ! add to the list of sets
    IF (.NOT.oldset) THEN
      CALL add_ele11 (elset, head, tail)
      nelms = nelms + 1 ! increased set counter for this element type
    END IF
    IF( oldset )DEALLOCATE(elset%stint)
    ALLOCATE(elset%stint((nstre+1)*ngaus,nelem))
    elset%stint = 0d0

  ELSE IF (TRIM(task) == 'RESTAR') THEN
    CALL new_ele11(elset)  !reserve memory for set
    ! read control parameters
    elset%sname = TRIM(sname)
    READ (51) nelem, nstre, ngaus, nbn, nreqs, narch
    READ (51) elset%shap, elset%stabs
    ALLOCATE(elset%stint((nstre+1)*ngaus,nelem))
    CALL rest11(nelem,nreqs,ntype,ngaus,nbn,elset%head,elset%tail,elset%ngrqs, &
                elset%nhead,elset%stint)
    elset%lside = .TRUE.
    CALL comm11 (0, nelem, nstre, nbn, ngaus, nreqs, narch, elsnam, elset)
    ! add to list of elements
    CALL add_ele11 (elset, head, tail)

  ELSE
    CALL runend('INPD11: NON-EXISTENT TASK .        ')
  ENDIF

  RETURN
  9999 CALL runen2('')

  END SUBROUTINE inpd11
