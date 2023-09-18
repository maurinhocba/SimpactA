 SUBROUTINE dump18( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

 TYPE (ele18_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs

 INTEGER (kind=4) :: i,j,size1,size2
 CHARACTER(len=midn):: straux
 TYPE (ele18), POINTER :: e

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'SOLAG'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, elset%gauss, elset%shell, &
       elset%locax, elset%bbar, elset%plstr, elset%angdf, elset%btscal, elset%small, elset%ngaus
 WRITE (50,ERR=9999) elset%gpc             !gauss-point coordinates
 e => elset%head                       !point to first element

 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT  !last element processed, EXIT loop
   ! dump variables associated to an element
   IF( ASSOCIATED(e%gausv) )THEN
     size1 = SIZE (e%gausv,1)
     size2 = SIZE (e%gausv,2)
   ELSE
     size1 = 0
     size2 = 0
   END IF
   WRITE (50,ERR=9999) e%numel,e%matno,e%lnods,size1,size2
   WRITE (50,ERR=9999) e%dvol,e%angle,e%cartd,e%stint
   IF( size1 > 0 )WRITE(50,ERR=9999) ((e%gausv(i,j),i=1,size1),j=1,size2 ) !Write internal variables
   IF( elset%shell ) WRITE(50,ERR=9999) e%nfdas,e%jacin                    !save values
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump18
