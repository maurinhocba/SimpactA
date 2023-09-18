 SUBROUTINE dump14( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************

 IMPLICIT NONE

 TYPE (ele14_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs
 INTEGER (kind=4) :: i,size1,size2,iel
 CHARACTER(len=midn):: straux
 TYPE (ele14), POINTER :: e

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'LBST'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, elset%logst, elset%lside, &
                     elset%gauss, elset%plstr, elset%angdf, elset%nonrg, elset%locax

 e => elset%head                       !point to first element
 iel = 0
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT  !last element processed, EXIT loop
   ! dump variables associated to an element
   iel = iel + 1
   IF( ASSOCIATED(e%gausv) )THEN
     size1 = SIZE (e%gausv,1)
     size2 = SIZE (e%gausv,2)
   ELSE
     size1 = 0
     size2 = 0
   END IF
   WRITE (50,ERR=9999) e%numel,e%matno,e%lnods,e%lside,size1,size2
   WRITE (50,ERR=9999) e%area1,e%lb,e%angle,e%a,e%b,e%ci,e%stra0,e%stra1,elset%stint(:,iel)
   IF( size1 > 0 )WRITE (50,ERR=9999) e%gausv
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump14
