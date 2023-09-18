 SUBROUTINE dump09(elset)
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

 TYPE(ele09_set) :: elset

 ! local variables
 TYPE(ele09), POINTER :: e
 CHARACTER(len=midn):: straux
 INTEGER (kind=4) :: ielem,i

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux
 straux = 'SHREV'
 WRITE (50,ERR=9999) straux,elset%sname
 WRITE (50,ERR=9999) elset%nelem
 WRITE (50,ERR=9999) elset%nnode, elset%ngaus, elset%nstre,   &
                     elset%axesc, elset%nreqs, elset%narch

 WRITE (50,ERR=9999) elset%posgp, elset%shape, elset%deriv, elset%weigh
 e => elset%head
 DO ielem=1,elset%nelem
   WRITE(50,ERR=9999) e%numel, e%matno, e%lnods, e%phil, e%jac, e%r0, e%stra0, &
                      e%stres, e%sede0, e%sedef
   IF( ASSOCIATED(e%ehist) )WRITE(50,ERR=9999) e%ehist, e%strap
   e => e%next
 END DO

 IF ( elset%nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,elset%nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump09
