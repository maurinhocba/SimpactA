 SUBROUTINE inpd19 (task, iwrit, elsnam, nelms)
 !   READ control DATA for element number 19 (TL PST++)
 USE ele19_db
 USE ctrl_db, ONLY : ntype
 IMPLICIT NONE

   ! dummy arguments
   CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
   CHARACTER(len=*),INTENT(IN):: task      ! requested task
   INTEGER (kind=4) :: nelms,   & ! number of element sets
                       iwrit      ! flag to echo data input
   ! local variables
   LOGICAL :: oldset
   INTEGER (kind=4) :: nreqs, narch, nelem, i

   TYPE (ele19_set), POINTER, SAVE  :: elset, anter


   IF (TRIM(task) == 'INPUT') THEN
     ! check if list of sets and set exists and initializes
     CALL srch_ele19 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm19 (1, nelem,  nreqs, narch, elsnam, elset)
       elset%lside = .FALSE.  !initializes flag to compute LSIDE
     ELSE                !set ELSNAM does not exist
       CALL new_ele19(elset)  !reserve memory for set
       CALL listen('INPD19')  !read a line
       elset%eulrf = exists('EULFRM') ! use a Eulerian Formulation
       IF (elset%eulrf .AND. ntype == 1) &
         CALL runen3('RESV19: EULER FORMULATION not allowed in plane stress')
       IF (exists('NOSWAP')) elset%swapc = .FALSE.  ! not swap alone elements in corners
       elset%angdf=getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
       nreqs=getint('NREQS ',0,' Gauss pt for stress time history..')
       !?      IF( nreqs > 0 )NULLIFY( elset%ngrqs )
       narch  =  0           !to check
       nelem  =  0           !new set, initializes number of elements
       IF(iwrit == 1)WRITE(lures,"(/,5X,'CONTROL PARAMETERS FOR 2-D TRIANGLE ' &
                     &       //,5X,'REQUIRED STRESS (NREQS) =',I10,/)",ERR=9999) nreqs
     END IF

     !  read new data or add to previous data
     CALL elmd19(nelem, elset%head, elset%tail, iwrit, elset%eulrf)
     elset%plstr = 0       ! do not compute plastic strains
     ALLOCATE(elset%eside(3,nelem)) !reserve memory for neighbor elements
     elset%eside = 0     !initializes
     IF (.NOT.oldset) CALL rdreqs (ngaus, nreqs, elset%ngrqs, iwrit )

     CALL comm19(0, nelem,  nreqs, narch, elsnam, elset)
     ! add to the list of sets
     IF (.NOT.oldset) THEN
       CALL add_ele19 (elset, head, tail)
       nelms = nelms + 1 ! increased set counter for this element type
     END IF

   ELSE IF (TRIM(task) == 'RESTAR') THEN

     CALL new_ele19(elset)  !reserve memory for set
     ! read control parameters
     elset%sname = elsnam
     READ (51) elset%nelem, elset%nreqs, elset%narch, &
               elset%gauss, elset%lside, elset%eulrf, elset%plstr, elset%angdf
     ALLOCATE(elset%eside(3,elset%nelem))
     READ(51) (elset%eside(:,i),i=1,elset%nelem)
     nreqs = elset%nreqs
     IF (nreqs > 0)THEN
       ALLOCATE(elset%ngrqs(nreqs))
       READ(51) (elset%ngrqs(i), i=1,nreqs)
     END IF
     ! restore list of elements
     CALL rest19 (elset%nelem,  elset%head, elset%tail)
     ! add to list of elements
     CALL add_ele19 (elset, head, tail)

   ELSE
     CALL runend('INPD19: NON-EXISTENT TASK .        ')
   END IF

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE inpd19
