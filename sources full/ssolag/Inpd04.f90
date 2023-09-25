 SUBROUTINE inpd04 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 04: SOLIDS solid-shell element

 USE ele04_db
 USE gvar_db,ONLY : fimpo,overw

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nreqs, narch, ngaus, nelem, i

 TYPE (ele04_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele04 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm04 (1, nelem,  nreqs, narch, elsnam, elset, ngaus)
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD04')  !read a line
     ngaus = getint('NGAUS ',2,' NUMBER OF GAUSS POINTS  (>= 2) ....')
     IF( ngaus < 2 )THEN
       WRITE(lures,"(' Minimumn number of gauss points is 2 ........')",ERR=9999)
       elset%ngaus = 2
     END IF
     elset%locax = 3       !default
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       IF( INT(param(i)) > 0 .AND. INT(param(i)) < 4 )THEN  !if a valid code
         elset%locax =  INT(param(i))
       ELSE  !keep default value and p
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") INT(param(i))
       END IF
     END IF
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%angdf = getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     elset%btscal   =getrea('BTSCAL',1d0,' Critical time increment scaler    ')
     elset%small = exists('SMALL')
     IF(elset%small) WRITE(lures,"(' Green strains will be used if possible')")
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements
     IF( exists('BETA  ',i)) THEN   !stabilization factors
       elset%beta = param(i:i+2)
       WRITE(lures,"(' Stabilization factors :',3f8.5)")elset%beta
     ELSE
       elset%beta = (/ 0.1d0, 0.05d0, 0.5d0 /)
     END IF
     IF( elset%beta(3) == 0d0 ) elset%beta(3) = 0.05
     !Initialize empty list Point both pointer no nothing
     CALL ini_ele04e (elset%head, elset%tail)
   END IF
   !  read new data or data to previous data
   CALL elmd04(nelem, elset%head, elset%tail, iwrit, ngaus)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele04 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   READ (51) nelem, nreqs, narch, ngaus,elset%gauss, elset%small,  &
             elset%plstr, elset%angdf, elset%locax, elset%btscal,  &
             elset%beta, elset%psg, elset%isg
   ! restore list of elements
   CALL rest04 (nelem, nreqs, elset%head, elset%tail, elset%ngrqs, ngaus )
   ! add to list of elements
   CALL add_ele04 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele04 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm04(1, nelem,  nreqs, narch, elsnam, elset, ngaus)
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     ALLOCATE (elset)       !reserve memory for set
     CALL ini_ele04e (elset%head, elset%tail)
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     READ (fimpo) nelem,ngaus,elset%angdf,elset%small,elset%locax,elset%beta
     nreqs = 0
     narch = 0
     elset%btscal   = 1d0   !Critical time increment scaler
   END IF
   ! restore list of elements
   CALL impo04 ( nelem, ngaus, elset%head, elset%tail)
   ! add to list of elements
   IF( .NOT.overw )THEN
     CALL add_ele04 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE
   CALL runend('INPD04: NON-EXISTENT TASK .        ')
 END IF
 CALL comm04(0, nelem,  nreqs, narch, elsnam, elset, ngaus)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd04
