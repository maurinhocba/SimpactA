 SUBROUTINE dump04( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

 TYPE (ele04_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs

 INTEGER (kind=4) :: i,size1
 CHARACTER(len=midn):: straux
 TYPE (ele04), POINTER :: e

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'PRISM'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, elset%gauss, &
       elset%plstr, elset%angdf, elset%btscal, elset%small, elset%nodvo
 e => elset%head                       !point to first element

 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT  !last element processed, EXIT loop
   ! dump variables associated to an element
   IF( ASSOCIATED(e%gausv) )THEN
     size1 = SIZE (e%gausv)
   ELSE
     size1 = 0
   END IF
   WRITE (50,ERR=9999) e%numel,e%matno,e%lnods,size1
   WRITE (50,ERR=9999) e%dvol,e%angle,e%cd,e%stint
   IF( size1 > 0 )WRITE(50,ERR=9999) (e%gausv(i),i=1,size1 ) !Write internal variables
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump04
