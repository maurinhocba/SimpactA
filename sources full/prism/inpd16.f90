 SUBROUTINE inpd16 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 16 (TL PRISM 6)

 USE ele16_db
 USE gvar_db,ONLY : fimpo,overw

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nelem, nreqs, narch, nn, ngaus, i

 TYPE (ele16_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele16 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm16 (1, nelem,  nreqs, narch, elsnam, elset, ngaus)

   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     CALL listen('INPD16')  !read a line
     nn = getint('NNODE ',6,' NUMBER OF ELEMENT NODES (6 ONLY)..')
     IF( nn /= nnode )CALL runend('PRISM: NNODE must be 6             ')
     ngaus = getint('NGAUS ',2,' NUMBER OF GAUSS POINTS  ( 2  ONLY)')
     IF( ngaus /= 2) CALL runend('PRISM: NGAUS must be 2             ')
     elset%bbar  = exists('BBAR ')
     IF(elset%bbar )WRITE(lures,"(' Average volumetric strain will be used')",ERR=9999)
     elset%quad = exists('QUAD')
     IF(elset%quad) WRITE(lures,"(9X,'QUADratic approximation will be used.')")
     elset%shell = exists('SHELL')
     IF(elset%shell)WRITE(lures,"(' Asumed transverse shear strains will be used')",ERR=9999)
     elset%lface = .FALSE.
     IF(elset%quad .OR. elset%shell ) THEN
       elset%locax = 3
     ELSE
       elset%locax = 0
     END IF
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       IF( elset%locax > 0 .AND. elset%locax < 4 )THEN  !if a valid code
         elset%locax =  INT(param(i))
       ELSE  !keep default value and p
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") INT(param(i))
       END IF
     END IF
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%angdf(1) =getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     elset%angdf(2) =getrea('BETA  ',0d0,' Second Euler Angle from X and Orth.')
     elset%angdf(3) =getrea('GAMMA ',0d0,' Third Euler Angle from X and Ortho')
     elset%btscal   =getrea('BTSCAL',1d0,' Critical time increment scaler    ')
     elset%small = exists('SMALL')
     IF(elset%small) WRITE(lures,"(' Green strains will be used if possible')")
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer no nothing
     CALL ini_ele16e (elset%head, elset%tail)
   END IF
   !  read new data or data to previous data
   CALL elmd16(nelem, elset%head, elset%tail, iwrit, ngaus, elset%quad, elset%shell)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele16 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   NULLIFY(elset%head)       !nullify head pointer
   ! read control parameters
   elset%sname = elsnam
   READ (51) elset%nelem, elset%btscal, elset%nreqs, elset%narch, elset%gauss, &
             elset%plstr, elset%angdf, elset%ngaus, elset%small, elset%shell,  &
             elset%quad, elset%locax, elset%bbar
   elset%lface = .TRUE.
   ! restore list of elements
   CALL rest16 (elset%nelem, elset%nreqs, elset%head, elset%tail, &
                elset%ngrqs, elset%ngaus, elset%quad, elset%shell  )
   ! add to list of elements
   CALL add_ele16 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele16 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm16(1, nelem,  nreqs, narch, elsnam, elset, ngaus)
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     ALLOCATE (elset)       !reserve memory for set
     CALL ini_ele16e (elset%head, elset%tail)
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     elset%lface = .FALSE.  !initializes flag to compute extended connectivities
     READ (fimpo) nelem,nn,ngaus,elset%angdf,elset%small,elset%quad, &
                  elset%bbar,elset%shell,elset%locax
     nreqs = 0
     narch = 0
     elset%btscal   = 1d0   !Critical time increment scaler
   END IF
   ! restore list of elements
   CALL impo16 ( nelem,ngaus, elset%head, elset%tail, elset%shell, elset%quad)
   ! add to list of elements
   IF( .NOT.overw )THEN
     CALL add_ele16 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE
   CALL runend('INPD16: NON-EXISTENT TASK .        ')
 END IF
 CALL comm16(0, nelem,  nreqs, narch, elsnam, elset, ngaus)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd16
