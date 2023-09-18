 SUBROUTINE inpd15 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 15 (TL BST++)

 USE ele15_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*), INTENT(IN):: elsnam     ! element set name
 CHARACTER(len=*), INTENT(IN):: task       ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets of this type
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset, logst
 INTEGER (kind=4) :: nreqs, nelem, narch
 CHARACTER(len=mnam) :: sname     ! element set name

 TYPE (ele15_set), POINTER, SAVE  :: elset, anter


 sname = elsnam
 IF (TRIM(task) == 'INPUT ') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele15 (head, anter, elset, sname, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm15 (1, nelem,  nreqs, narch, sname, elset, logst)
     elset%lside = .FALSE.  !initializes flag to compute LSIDE
   ELSE                !set ELSNAM does not exist
     CALL new_ele15(elset)
     CALL listen('INPD15')  !read a line
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%angdf=getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     logst = .NOT.exists('SMALL ')
     IF( exists('GSHEAR ')) elset%shear = 1
     IF( exists('NSHEAR ')) elset%shear =-1
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     IF(iwrit == 1)WRITE(lures,"(/,5X,'CONTROL PARAMETERS FOR SHELL ELEMENT' &
                   &       //,5X,'REQUIRED STRESS (NREQS) =',I10,  &
                   &       //,5X,'USE G-L strains (SMALL) =',L10 ,/)",ERR=9999)&
                   nreqs,.NOT.logst

     !Initialize empty list Point both pointer to nothing
     CALL ini_ele15e (elset%head, elset%tail)
   END IF
   !  read new data or add to previous data
   CALL elmd15(nelem, elset%head, elset%tail, iwrit)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL comm15(0, nelem,  nreqs, narch, sname, elset, logst)

   CALL elmd15r(elset%nrf, elset%rhead, elset%rtail, iwrit)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele15 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   CALL new_ele15(elset)      !initializes a list
   ! read control parameters
   elset%sname = sname
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%nrf,  elset%logst, &
             elset%lside, elset%gauss, elset%plstr, elset%angdf , elset%shear
   ! restore list of elements
   ALLOCATE (elset%stint(10,elset%nelem))         !initializes a list
   CALL rest15 (elset%nelem,  elset%nreqs, elset%head,  elset%tail, &
                elset%ngrqs,  elset%nrf,   elset%rhead, elset%rtail, elset%stint )
   ! add to list of elements
   CALL add_ele15 (elset, head, tail)

 ELSE
   CALL runend('INPD15: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd15
