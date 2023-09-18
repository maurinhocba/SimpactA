 SUBROUTINE inpd20 (task, iwrit, elsnam, nelms)
 !   READ control DATA for element number 20 (TL PST++)
 USE ele20_db
 USE ctrl_db, ONLY : ntype
 USE gvar_db,ONLY : maxsdv,fimpo,overw
 IMPLICIT NONE

   ! dummy arguments
   CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
   CHARACTER(len=*),INTENT(IN):: task      ! requested task
   INTEGER (kind=4) :: nelms,   & ! number of element sets
                       iwrit      ! flag to echo data input
   ! local variables
   LOGICAL :: oldset
   INTEGER (kind=4) :: nreqs, narch, nelem, i
   CHARACTER(len=mnam) :: sname     ! element set name

   TYPE (ele20_set), POINTER, SAVE  :: elset, anter

   sname = elsnam
   IF (TRIM(task) == 'INPUT') THEN
     ! check if list of sets and set exists and initializes
     CALL srch_ele20 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm20 (1, nelem,  nreqs, narch, elsnam, elset)
       elset%lside = .FALSE.  !initializes flag to compute LSIDE
     ELSE                !set ELSNAM does not exist
       CALL new_ele20(elset)  !reserve memory for set
       CALL listen('INPD20')  !read a line
       elset%eulrf = exists('EULFRM') ! use a Eulerian Formulation
       IF (elset%eulrf .AND. ntype == 1) &
         CALL runen3('RESV20: EULER FORMULATION not allowed in plane stress')
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
     CALL elmd20(nelem, elset%head, elset%tail, iwrit, elset%eulrf)
     elset%plstr = 0       ! do not compute plastic strains
     ALLOCATE(elset%eside(3,nelem)) !reserve memory for neighbor elements
     elset%eside = 0     !initializes
     IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

     CALL comm20(0, nelem,  nreqs, narch, elsnam, elset)
     ! add to the list of sets
     IF (.NOT.oldset) THEN
       CALL add_ele20 (elset, head, tail)
       nelms = nelms + 1 ! increased set counter for this element type
     END IF

   ELSE IF (TRIM(task) == 'RESTAR') THEN

     CALL new_ele20(elset)  !reserve memory for set
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
     CALL rest20 (elset%nelem,  elset%head, elset%tail)
     ! add to list of elements
     CALL add_ele20 (elset, head, tail)

   ELSE IF (TRIM(task) == 'IMPORT') THEN
     ! check if list of sets and set exists and initializes
     IF( overw )THEN
       CALL srch_ele20 (head, anter, elset, elsnam, oldset)
       IF (oldset) THEN    !if set exists
         CALL comm20 (1, nelem,  nreqs, narch, sname, elset)
         !elset%lside = .FALSE.   !initializes flag to compute LSIDE
       ELSE
         CALL runen2(' Old set to overwrite does not exist')
       END IF
       READ (fimpo)
     ELSE
       CALL new_ele20(elset)
       elset%sname = sname
       READ (fimpo) nelem,elset%angdf,elset%eulrf
       nreqs = 0
       narch = 0
     END IF
     ! restore list of elements
     CALL impo20 ( nelem,elset%head, elset%tail)
     elset%origl = .FALSE. ! node label are changed to internal numeration
     CALL comm20 (0, nelem,  nreqs, narch, sname, elset)
     ! add to list of elements
     IF( .NOT.overw )THEN
       CALL add_ele20 (elset, head, tail)
       nelms = nelms + 1 ! increased set counter for this element type
     END IF

   ELSE
     CALL runend('INPD20: NON-EXISTENT TASK .        ')
   END IF

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE inpd20
