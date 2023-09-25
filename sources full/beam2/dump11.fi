 SUBROUTINE dump11(elset)
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 USE ctrl_db, ONLY : ntype
 IMPLICIT NONE

 TYPE(ele11_set) :: elset

 ! local variables
 TYPE(ele11), POINTER :: e
 TYPE(nodei), POINTER :: ni
 INTEGER (kind=4) :: ielem,i,nvar,nlay
 CHARACTER(len=midn):: straux
 LOGICAL :: fcur

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'BEAM2'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem
 WRITE (50,ERR=9999) elset%nelem, elset%nstre, elset%ngaus, elset%nbn, elset%nreqs, &
                     elset%narch
 WRITE (50,ERR=9999) elset%shap,elset%stabs

 e => elset%head
 DO ielem=1,elset%nelem
   fcur = ASSOCIATED( e%cur )
   WRITE(50,ERR=9999)e%numel, e%matno, e%lnods, e%l0, e%cab, e%lambd, elset%stint(:,ielem), fcur
   IF(ntype == 3 )WRITE(50,ERR=9999)e%ro
   IF( ASSOCIATED(e%gausv) )THEN
     nvar = SIZE(e%gausv,1)
     nlay = SIZE(e%gausv,2)
     WRITE(50,ERR=9999)nvar,nlay
     WRITE(50,ERR=9999)e%gausv
   ELSE
     nvar = 0
     nlay = 0
     WRITE(50,ERR=9999)nvar,nlay
   END IF
   IF( fcur )WRITE(50,ERR=9999)e%cur

   e => e%next
 END DO

 IF ( elset%nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,elset%nreqs)

 IF ( elset%nbn > 0 )THEN
   ni => elset%nhead
   DO i=1,elset%nbn
     WRITE (50,ERR=9999) ni%nn,ni%ba
     WRITE (50,ERR=9999) ni%lnods,ni%alph0
     ni => ni%next
   END DO
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump11
