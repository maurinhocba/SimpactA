 SUBROUTINE inpd13 (task, iwrit, elsnam, nelms)

 !   READ control DATA for element number 13 (TL BST++)

 USE ele13_db
 USE gvar_db,ONLY : maxsdv,fimpo,overw

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

 TYPE (ele13_set), POINTER, SAVE  :: elset, anter


 sname = elsnam
 IF (TRIM(task) == 'INPUT ') THEN
   ! check if list of sets and set exists and initializes
   CALL srch_ele13 (head, anter, elset, sname, oldset)
   IF (oldset) THEN    !if set exists
     CALL comm13 (1, nelem,  nreqs, narch, sname, elset, logst)
     elset%lside = .FALSE.  !initializes flag to compute LSIDE
   ELSE                !set ELSNAM does not exist
     CALL new_ele13(elset)
     CALL listen('INPD13')  !read a line
     nreqs=getint('NREQS ',0,' GAUSS PT FOR STRESS TIME HISTORY..')
     IF( nreqs > 0 )NULLIFY( elset%ngrqs )
     elset%angdf=getrea('ANGLE ',0d0,' Default angle between X1 and Ort_1')
     IF( exists('LOCAL  ',i)) THEN   !local axis definition
       elset%locax =  INT(param(i))
       IF( elset%locax < 1 .OR. elset%locax > 3 )THEN
         WRITE(lures,"(/,5X,'Error in the definition of local system', / &
                         5X,'Invalid axis: ',i3,' Default value used (3)')") elset%locax
         elset%locax = 3
       END IF
     END IF
     logst = .NOT.exists('SMALL ')
     IF( exists('GSHEAR ')) elset%shear = 1
     IF( exists('NSHEAR ')) elset%shear =-1
     narch  =  0           !to check
     nelem  =  0           !new set, initializes number of elements

     IF(iwrit == 1)WRITE(lures,"(/,5X,'CONTROL PARAMETERS FOR SHELL ELEMENT' &
                   &       //,5X,'REQUIRED STRESS (NREQS) =',I10,  &
                   &       //,5X,'USE G-L strains (SMALL) =',L10,  &
                   &       //,5X,'Include Shear st(SHEAR) =',i10,/)",ERR=9999)&
                   nreqs,.NOT.logst,elset%shear

     !Initialize empty list Point both pointer to nothing
     CALL ini_ele13e (elset%head, elset%tail)
   END IF
   !  read new data or add to previous data
   CALL elmd13(nelem, elset%head, elset%tail, iwrit, elset%shear)
   elset%plstr = 0     ! do not compute plastic strains
   IF (.NOT.oldset) CALL rdreqs ( 1 ,nreqs, elset%ngrqs, iwrit )

   CALL comm13(0, nelem,  nreqs, narch, sname, elset, logst)
   ! add to the list of sets
   IF (.NOT.oldset) THEN
     CALL add_ele13 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF
   IF( ASSOCIATED(elset%stint) )DEALLOCATE(elset%stint)
   ALLOCATE(elset%stint(10,nelem))
   elset%stint = 0d0

 ELSE IF (TRIM(task) == 'RESTAR') THEN

   CALL new_ele13(elset)
   ! read control parameters
   elset%sname = sname
   READ (51) elset%nelem, elset%nreqs, elset%narch, elset%nbs,  elset%logst, &
             elset%lside, elset%gauss, elset%plstr, elset%angdf, elset%shear, elset%locax
   ! restore list of elements
   ALLOCATE (elset%stint(10,elset%nelem))         !(10) initializes a list
   CALL rest13 (elset%nelem,  elset%nreqs, elset%head,  elset%tail, elset%shear, &
                elset%ngrqs,  elset%nbs,   elset%bhead, elset%btail, elset%stint,&
                elset%factors,elset%ninv, elset%moments  )
   ! add to list of elements
   CALL add_ele13 (elset, head, tail)

 ELSE IF (TRIM(task) == 'IMPORT') THEN
   ! check if list of sets and set exists and initializes
   IF( overw )THEN
     CALL srch_ele13 (head, anter, elset, elsnam, oldset)
     IF (oldset) THEN    !if set exists
       CALL comm13 (1, nelem,  nreqs, narch, sname, elset, logst)
       !elset%lside = .FALSE.   !initializes flag to compute LSIDE
     ELSE
       CALL runen2(' Old set to overwrite does not exist')
     END IF
     READ (fimpo)
   ELSE
     CALL new_ele13(elset)
     elset%sname = sname
     READ (fimpo) nelem,elset%angdf,logst,elset%locax
     nreqs = 0
     narch = 0
     ALLOCATE(elset%stint(10,nelem))
     elset%stint = 0d0
   END IF
   ! restore list of elements
   CALL impo13 ( nelem,elset%head, elset%tail)
   elset%origl = .FALSE. ! node label are changed to internal numeration
   CALL comm13 (0, nelem,  nreqs, narch, sname, elset, logst)
   ! add to list of elements
   IF( .NOT.overw )THEN
     CALL add_ele13 (elset, head, tail)
     nelms = nelms + 1 ! increased set counter for this element type
   END IF

 ELSE
   CALL runend('INPD13: NON-EXISTENT TASK .        ')
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpd13
