 SUBROUTINE immtsc (nf,nm,ns,mn,sn)
 !
 !   import materials and sections
 !
 USE mat_dba

 IMPLICIT NONE

 INTEGER(kind=4) :: nf,nm,ns
 INTEGER(kind=4) :: mn(2,nm),sn(2,ns)

 LOGICAL :: found,equal,exist(0:99)  !flags
 INTEGER (kind=4) :: i,j,k
 TYPE(mater), POINTER :: mat,matn
 TYPE(section), POINTER :: sec,secn
 TYPE(curve), POINTER :: cur
 TYPE(postv), POINTER :: postp,post1
 REAL (kind=8) :: tol
 TYPE(flc_tp),POINTER:: flc

 !   ************** READ   M A T E R I A L S  ****************************

 !  first task, see material labels
 exist = .FALSE.
 mat => mhead
 DO i=1,numats
   j = mat%matno
   IF( j < 100 )exist(j) = .TRUE.
   mat => mat%next
 END DO

 DO i=1,nm  !for each material stored

   CALL ini_mat( matn )   !get memory for a material
   READ (nf) matn%matno, matn%mtype, matn%matdef !read all the control variables
   ALLOCATE( matn%prope(matn%matdef(9)), &      !get memory
             matn%propp(matn%matdef(10)), &
             matn%props(matn%matdef(11)))
   READ (nf) matn%prope, matn%propp, matn%props  !read all the properties
   ! read associated curves if exists
   IF( matn%matdef(12) > 0 )THEN
     DO k=1,matn%matdef(12)   !for each curve
       ALLOCATE(cur)
       READ(nf)cur%np                !number of points in the curve
       ALLOCATE( cur%val(3,cur%np) ) !get memory
       READ(nf)cur%val               !curve values
       CALL add_cur(cur,matn%chead,matn%ctail) !add curve to the list
     END DO
   END IF
   ! compare material read with materials in data base
   equal = .FALSE.
   mat => mhead
   DO j=1,numats
     equal = matn%mtype == mat%mtype
     DO k=1,12
       equal = equal .AND. matn%matdef(k) == mat%matdef(k)
     END DO
     IF(.NOT. equal ) CYCLE
     DO k=1,matn%matdef(9)
       tol = MAX(0.0000001, ABS(matn%prope(k))/10000000)
       equal = equal .AND. ABS(matn%prope(k) - mat%prope(k)) < tol
     END DO
     IF(.NOT. equal ) CYCLE
     DO k=1,matn%matdef(10)
       tol = MAX(0.0000001, ABS(matn%propp(k))/10000000)
       equal = equal .AND. ABS(matn%propp(k) - mat%propp(k)) < tol
     END DO
     IF(.NOT. equal ) CYCLE
     DO k=1,matn%matdef(11)
       tol = MAX(0.0000001, ABS(matn%props(k))/10000000)
       equal = equal .AND. ABS(matn%props(k) - mat%props(k)) < tol
     END DO
     ! curves not compared yet  !!
     IF( equal ) EXIT           !identical material found
     mat => mat%next            !next material in database
   END DO
   mn(1:2,i) = matn%matno !remember material name
   IF( equal )THEN        !if an equal material already exists
     mn(2,i) = mat%matno  !keep new name for this material
     ! release memory of the material
     DEALLOCATE( matn%prope,matn%propp,matn%props )
     IF( matn%matdef(12) > 0 )THEN   !
       cur => matn%chead
       DO k=1,matn%matdef(12)   !for each curve
         DEALLOCATE(cur%val)
         cur => cur%next
       END DO
       NULLIFY(matn%chead)
       NULLIFY(matn%ctail)
     END IF
     DEALLOCATE (matn)
   ELSE
     k = matn%matno       !material label
     IF( k < 100 )THEN   ! for low name materials
       IF( exist(matn%matno) )THEN
         DO j=0,99         !find the first not used name
           IF( exist(j) )CYCLE
           exist(j) = .TRUE.
           matn%matno = j   !assign new name
           mn(2,i) = j      !keep new name for this material
           EXIT             !exit loop
         END DO
       END IF
     ELSE                 ! for high name materials
       DO j=50,99         !find the first not used name
         IF( exist(j) )CYCLE
         exist(j) = .TRUE.
         matn%matno = j   !assign new name
         mn(2,i) = j      !keep new name for this material
         EXIT             !exit loop
       END DO
     END IF
     CALL add_mat( matn, mhead, mtail )  !add material to the list
     numats = numats + 1
   END IF
 END DO

 IF( nm > 0 )THEN
   IF( ASSOCIATED(pmats) )DEALLOCATE(pmats)
   ALLOCATE (pmats(numats))
   mat => mhead
   DO i=1,numats
     pmats(i)%p => mat
     mat => mat%next
   END DO
 END IF


 !   ************** READ   S E C T I O N S  ****************************

 !  first task, see material labels
 exist = .FALSE.
 sec => shead
 DO i=1,nusect
   j = sec%secno
   IF( j < 100 )exist(j) = .TRUE.
   sec => sec%next
 END DO

 IF( nflc == 0 ) CALL ini_flc_db( )

 DO i=1,ns  !for each section stored

   ALLOCATE( secn )   !get memory for a section
   READ (nf) secn%secno, secn%secty, secn%mabas, secn%secdef  !read all the control variables

   ALLOCATE( secn%iprop(secn%secdef(1)), &      !get memory
             secn%rprop(secn%secdef(2)))

   READ (nf) secn%iprop, secn%rprop  !read all the properties

   ! read user-defined postprocess information, if exists
   IF( secn%secdef(4) > 0)THEN
     ALLOCATE (secn%postp)               !get memory for first variable
     postp => secn%postp
     k = 0                               !initializes counter
     DO                                  !loop
       k = k+1                           !increase loop counter
       READ (nf) postp%type,postp%dim,postp%name   !read data
       IF( k == secn%secdef(4) )EXIT     !if all variables read, exit
       ALLOCATE (post1)                  !get memory for next variable
       postp%next => post1               !put in the list
       postp => postp%next               !go to next
     END DO
   END IF

   IF( secn%secdef(5) >= 3)THEN  !if the section may have FLC
     IF( secn%iprop(3) > 0)THEN  !If the FLC label is greater than 0
       CALL new_flc(flc)
       READ(nf) flc%lbl,flc%sttyp,flc%npt                !Number of points of the FLC curve
       READ(nf) flc%LmM, flc%LmPS, flc%LmLS, flc%LmT     !read  strains limits
       ALLOCATE(flc%cv(2,flc%npt))
       READ(nf) (flc%cv(1:2,j),j=1,flc%npt)              !read  strains coordinates of the FLC
       CALL add_flc(flc,hpflc,tpflc)
     END IF
   END IF
   ! find material used
   j = 0
   DO
    j = j+1
    IF( secn%mabas == mn(1,i) )EXIT
   END DO
   secn%mabas = mn(2,i)    !assign new name

   ! compare section read with sections in data base
   equal = .FALSE.
   sec => shead
   DO j=1,nusect
     equal = secn%secty == sec%secty  .AND. secn%mabas == sec%mabas
     IF(.NOT. equal ) CYCLE
     DO k=1,5
       equal = equal .AND. secn%secdef(k) == sec%secdef(k)
     END DO
     IF(.NOT. equal ) CYCLE
     DO k=1,secn%secdef(1)
       equal = equal .AND. (secn%iprop(k) == sec%iprop(k))
     END DO
     IF(.NOT. equal ) CYCLE
     DO k=1,secn%secdef(2)
       tol = MAX(0.0000001, ABS(secn%rprop(k))/10000000)
       equal = equal .AND. ABS(secn%rprop(k) - sec%rprop(k)) < tol
     END DO
     IF(.NOT. equal ) CYCLE
     ! post-process information not compared yet !!
     IF( equal ) EXIT           !identical section found
     sec => sec%next            !next section in database
   END DO
   sn(1:2,i) = secn%secno !remember section name
   IF( equal )THEN        !if an equal section already exists
     sn(2,i) = sec%secno  !keep new name for this section
     ! release memory of the section
     DEALLOCATE( secn%iprop,secn%rprop )
     IF( secn%secdef(4) > 0 ) NULLIFY(secn%postp)
     DEALLOCATE (secn)
   ELSE
     k = secn%secno       !section label
     IF( k < 100 )THEN   !only for low name sections
       IF( exist(secn%secno) )THEN
         DO j=0,99         !find the first not used name
           IF( exist(j) )CYCLE
           exist(j) = .TRUE.
           secn%secno = j   !assign new name
           sn(2,i) = j      !keep new name for this section
           EXIT             !exit loop
         END DO
       END IF
     ELSE                ! for high name sections
       DO j=50,99         !find the first not used name
         IF( exist(j) )CYCLE
         exist(j) = .TRUE.
         secn%secno = j   !assign new name
         sn(2,i) = j      !keep new name for this section
         EXIT             !exit loop
       END DO
     END IF
     CALL mat_search (secn%mabas,found,secn%mtbas)
     CALL add_sect( secn, shead, stail )  !add section to the list
     nusect = nusect + 1
   END IF
 END DO
 IF( ns > 0 )THEN
   IF (ASSOCIATED(psecs)) DEALLOCATE (psecs)
   ALLOCATE (psecs(nusect))
   sec => shead
   DO i=1,nusect
     psecs(i)%p => sec
     sec => sec%next
   END DO
 END IF

 RETURN
 END SUBROUTINE immtsc
