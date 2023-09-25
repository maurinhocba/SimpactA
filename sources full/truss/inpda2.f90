 SUBROUTINE inpda2(task,nelem,iwrit,elsnam,nelms)
 !******************************************************************
 !
 !*** READ control DATA for 2-node 2/3-d truss element
 !
 !******************************************************************

 USE ele02_db

 IMPLICIT NONE

 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: iwrit,   & ! flag to echo data input
&                    nelem,   & ! number of elements
&                    nelms      ! number of element sets of this type

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nreqs,narch

 TYPE (ele02_set), POINTER :: elset,anter


 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele02 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL commv2 (1,nelem,nreqs,narch,elsnam,elset)
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     CALL listen('INPDA2')
     nreqs = getint('NREQS ',0 ,'!Gauss pt for stress time history..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     !Initialize empty list Point both pointer to nothing
     CALL ini_ele02e (elset%head, elset%tail)
   END IF
   CALL elmda2(nelem,elset%head,elset%tail,iwrit)

   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL commv2 (0,nelem,nreqs,narch,elsnam,elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele02 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)       !reserve memory for set
   CALL ini_ele02e (elset%head, elset%tail)
   READ (51) nelem,nreqs,narch
   CALL resta2 (nelem,nreqs,elset%head, elset%tail,elset%ngrqs)
   CALL add_ele02 (elset, head, tail)

 ELSE
   CALL runend('INPDA2: NON-EXISTENT TASK .        ')
 END IF

 CALL commv2 (0,nelem,nreqs,narch,elsnam,elset)

 RETURN
 END SUBROUTINE inpda2
