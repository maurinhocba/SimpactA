 SUBROUTINE dumpi1 (elset)
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE
 TYPE(ele01_set), POINTER :: elset          !element set

 INTEGER (kind=4) :: i,ielem,nvare
 CHARACTER(len=midn):: straux
 TYPE(ele01), POINTER :: e          !element set

 nvare = SIZE(elset%head%gausv)     !

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux
 straux = 'SPOT'
 WRITE (50,ERR=9999) straux,elset%sname
 WRITE (50,ERR=9999) elset%nelem
 WRITE (50,ERR=9999) elset%nnode,elset%nelem, elset%nreqs, elset%narch, nvare

 e => elset%head
 DO ielem=1,elset%nelem
   WRITE (50,ERR=9999) e%lnods , e%numel, e%matno, e%rot
   WRITE (50,ERR=9999) e%gausv
   e => e%next
 END DO
 IF ( elset%nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,elset%nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dumpi1
