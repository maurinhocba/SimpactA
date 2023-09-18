 SUBROUTINE inpd04 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 04 (TL CST++)

 USE ele04_db

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nelem, nreqs, narch, nn

 TYPE (ele04_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT ') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele04 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm04 (1, nelem,  nreqs, narch, elsnam, elset)
     elset%lside = .NOT.elset%nodvo   !initializes flag to compute LSIDE
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD04')  !read a line
     nn = getint('NNODE ',4,' NUMBER OF ELEMENT NODES (4 ONLY)..')
     IF( nn /= nnode )CALL runend('TETRA: NNODE must be 4             ')
     nreqs=getint('NREQS ',0,' Gauss pt for stress time history..')
     elset%angdf(1) =getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     elset%angdf(2) =getrea('BETA  ',0d0,' Second Euler Angle from X and Orth')
     elset%angdf(3) =getrea('GAMMA ',0d0,' Third Euler Angle from X and Ortho')
     elset%nodvo  = .NOT.exists('STANDA')
     IF( elset%nodvo  )THEN
       WRITE(lures,"(/,5X,'Volume Approximation based on neighbours will be used')")
       nodvo = .TRUE.
       stab  =getrea('STABIL',stab,' Stabiliztion factor               ')
     ELSE
       WRITE(lures,"(/,5X,'STANDARD Approximation Strain) will be used')")
     END IF
     elset%btscal = getrea('BTSCAL',1d0,' Critical time increment scaler    ')
     elset%lside  = .NOT.elset%nodvo   !initializes flag to compute LSIDE
     IF(exists('SMALL'))THEN
       elset%small = .TRUE.
       WRITE(lures,"(' Green strains will be used if possible')")
     ELSE
       elset%small = .FALSE.
     END IF

     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer to nothing
     CALL ini_ele04e (elset%head, elset%tail)
   END IF
   !  read new data or add to previous data
   CALL elmd04( nelem, elset%head, elset%tail, iwrit )
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   CALL comm04(0, nelem,  nreqs, narch, elsnam, elset)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele04 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   !  read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%gauss, &
       elset%plstr, elset%angdf, elset%btscal, elset%small, elset%nodvo
   CALL rest04 (elset%nelem,  elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs)
   ! add to list of elements
   CALL add_ele04 (elset, head, tail)

 ELSE
   CALL runend('INPD04: NON-EXISTENT TASK .        ')
 END IF
 RETURN

 END SUBROUTINE inpd04
