 SUBROUTINE inpd25 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 25 (TL BSQ)

 USE ele25_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*), INTENT(IN):: elsnam     ! element set name
 CHARACTER(len=*), INTENT(IN):: task       ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset, logst
 INTEGER (kind=4) :: nreqs, narch, nelem, i
 CHARACTER(len=mnam) :: sname     ! element set name

 TYPE (ele25_set), POINTER, SAVE :: elset, anter

 sname = elsnam
 IF (TRIM(task) == 'INPUT ') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele25 (head, anter, elset, sname, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm25 (1, nelem,  nreqs, narch, sname, elset, logst)
     elset%lside = .FALSE.  !initializes flag to compute LSIDE
   ELSE                !set ELSNAM does not exist
     CALL new_ele25(elset)
     CALL listen('INPD25')  !read a line
     nreqs = getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     elset%angdf = getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 3
       END IF
     END IF
     elset%stabs = getrea('STABS ',0.15D0,' Membrane Stabilization coeff. ... ')
     elset%stabb = getrea('STABB ',0.05D0,' Bending Stabilization coeff. .... ')
     logst = .NOT.exists('SMALL ')
     IF (logst) WRITE(lures,"(T10,'Large Strain formulation (Hencky) Will be used')",ERR=9999)
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     IF(iwrit == 1)WRITE(lures,"(/,5X,'CONTROL PARAMETERS FOR SHELL ELEMENT' &
                   &       //,5X,'REQUIRED STRESS (NREQS) =',I10,/)",ERR=9999) nreqs
   END IF
   !  read new data or add to previous data
   CALL elmd25(nelem, elset%head, elset%tail, iwrit)
   IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
   ALLOCATE(elset%stint(10,nelem))
   elset%stint = 0d0
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL comm25(0, nelem,  nreqs, narch, sname, elset, logst)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele25 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = sname
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%nbs, elset%logst, &
             elset%lside, elset%gauss, elset%plstr, elset%angdf, elset%stabs, elset%stabb
   ! restore list of elements
   ALLOCATE (elset%stint(10,elset%nelem))         !initializes a list
   CALL rest25 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%nbs, elset%bhead, elset%btail, elset%stint )
   ! add to list of elements
   CALL add_ele25 (elset, head, tail)

 ELSE
   CALL runend('INPD25: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')

 END SUBROUTINE inpd25
