 SUBROUTINE dumpi2 (elset)
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

 TYPE(ele02_set), POINTER :: elset          !element set

 INTEGER (kind=4) :: i,ielem
 CHARACTER(len=midn):: straux
 TYPE(ele02), POINTER :: e          !element set

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux
 straux = 'TRUSS'
 WRITE (50,ERR=9999) straux,elset%sname
 WRITE (50,ERR=9999) elset%nelem
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch

 e => elset%head
 DO ielem=1,elset%nelem
   WRITE (50,ERR=9999) e%lnods , e%numel, e%matno
   WRITE (50,ERR=9999) e%gausv
   e => e%next
 END DO
 IF ( elset%nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,elset%nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dumpi2
