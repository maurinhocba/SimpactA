 SUBROUTINE inpd18 (task, nel, iwrit, elsnam, nelms)

 !   READ control DATA for element number 18 (TL QUAD 4)

 USE ele18_db
 USE gvar_db,ONLY : maxsdv,fimpo,overw

 IMPLICIT NONE
 ! dummy arguments
 CHARACTER(len=*),INTENT(IN):: elsnam    ! element set name
 CHARACTER(len=*),INTENT(IN):: task      ! requested task
 INTEGER (kind=4) :: nelms,   & ! number of element sets
                     nel,     & ! number of elements in the set
                     iwrit      ! flag to echo data input

 ! local variables
 LOGICAL :: oldset
 INTEGER (kind=4) :: nelem, nreqs, narch, nn, ngaus, i

 TYPE (ele18_set), POINTER, SAVE  :: elset, anter

 IF (TRIM(task) == 'INPUT') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele18 (head, anter, elset, elsnam, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm18 (1, nelem,  nreqs, narch, elsnam, elset, ngaus)
     nel = nel - nelem
   ELSE                !set ELSNAM does not exist
     ALLOCATE (elset)       !reserve memory for set
     !initializes flags
     elset%gauss = .FALSE.  !flag to compute Gauss constants
     elset%shell = .FALSE.
     elset%check = .FALSE.
     elset%bbar  = .TRUE.
     elset%locax = 0
     CALL listen('INPD18')  !read a line
     nn = getint('NNODE ',8,' NUMBER OF ELEMENT NODES (8 ONLY)..')
     IF( nn /= nnode )CALL runend('SOLAG: NNODE must be 8             ')
     ngaus = getint('NGAUS ',8,' NUMBER OF GAUSS POINTS  (4/8 ONLY)')
     IF( ngaus == 2 ) ngaus = 8
     IF( ngaus /= 4 .AND. ngaus /= 8 )CALL runend (' ERROE IN NUMBER OF GAUSS POINTS  (4/8 ONLY)')
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF(exists('SHELL'))THEN  !shell version
       elset%shell = .TRUE.
       elset%locax = 3
       WRITE(lures,"(9x,'SHELL option: Assumed transverse shear strains will be used')",ERR=9999)
       IF(exists('CHECK'))THEN  !CHECK option is activated
         elset%check = .TRUE.
         WRITE(lures,"(9x,'CHECK option: Connectivities order will be checked')",ERR=9999)
       END IF
     ELSE
       elset%locax = 0        !global XYZ system
     END IF
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 0
       END IF
     END IF
     elset%angdf(1) =getrea('ALPHA ',0d0,' First Euler Angle from X and Ortho')
     elset%angdf(2) =getrea('BETA  ',0d0,' Second Euler Angle from X and Orth.')
     elset%angdf(3) =getrea('GAMMA ',0d0,' Third Euler Angle from X and Ortho')
     elset%btscal   =getrea('BTSCAL',1d0,' Critical time increment scaler    ')
     IF(exists('SMALL'))THEN
       elset%small = .TRUE.
       WRITE(lures,"(' Green strains will be used if possible')",ERR=9999)
     ELSE
       elset%small = .FALSE.
     END IF
     IF(exists('NOBBAR'))THEN
       elset%bbar  = .FALSE.
       WRITE(lures,"(' BBAR formulation dissabled ')",ERR=9999)
     END IF

     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     !Initialize empty list Point both pointer no nothing
     CALL ini_ele18e (elset%head, elset%tail)
     ALLOCATE( elset%gpc(3,ngaus))
   END IF
   !  read new data or data to previous data
   CALL elmd18(nelem, elset%head, elset%tail, iwrit, ngaus,elset%shell)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs (ngaus,nreqs, elset%ngrqs, iwrit )

   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele18 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   nel = nel + nelem

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   ALLOCATE (elset)          !initializes a list
   elset%check = .FALSE.
   CALL ini_ele18e (elset%head, elset%tail)     !nullify head pointer
   ! read control parameters
   READ (51) nelem, nreqs, narch, elset%gauss, elset%shell, elset%locax, elset%bbar, &
             elset%plstr, elset%angdf, elset%btscal, elset%small, ngaus
   ALLOCATE( elset%gpc(3,ngaus))
   READ (51) elset%gpc
   ! restore list of elements
   CALL rest18 (nelem, nreqs, elset%head, elset%tail, &
                elset%ngrqs, ngaus, elset%shell  )
   ! add to list of elements
   CALL add_ele18 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele18 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm18(1, nelem,  nreqs, narch, elsnam, elset, ngaus)
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     ALLOCATE (elset)       !reserve memory for set
     CALL ini_ele18e (elset%head, elset%tail)
     elset%gauss = .FALSE.  !initializes flag to compute Gauss constants
     READ (fimpo) nelem,ngaus,elset%angdf,elset%small,elset%shell,elset%bbar,elset%locax
     nreqs = 0
     narch = 0
     ALLOCATE( elset%gpc(3,ngaus))
     elset%btscal   = 1d0   !Critical time increment scaler
     elset%check = .FALSE.  !do not check connectivities
   END IF
   ! restore list of elements
   CALL impo18 ( nelem,ngaus, elset%head, elset%tail, elset%shell)
   elset%check = .FALSE.
   ! add to list of elements
   IF( .NOT.overw )THEN
     CALL add_ele18 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE
   CALL runend('INPD18: NON-EXISTENT TASK .        ')
 END IF
 CALL comm18(0, nelem,  nreqs, narch, elsnam, elset, ngaus)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd18
