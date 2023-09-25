 SUBROUTINE dump13( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE
 TYPE (ele13_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs

 INTEGER (kind=4) :: i,size1,size2,iel
 TYPE (ele13), POINTER :: e
 TYPE (sideb), POINTER :: sb
 CHARACTER(len=midn):: straux

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'NBST'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, elset%nbs,  elset%logst, &
                     elset%lside, elset%gauss, elset%plstr, elset%angdf, elset%shear, elset%locax

 ! dump branching sides data base
 IF( elset%nbs > 0 )THEN
   sb => elset%bhead
   DO i=1,elset%nbs
     ! dump variables associated to the side
     WRITE (50,ERR=9999) sb%nn
     WRITE (50,ERR=9999) sb%lnods,sb%alph0,sb%fc,sb%c
     !WRITE (50,ERR=9999) sb%gamma,sb%bb  !unnecessary
     sb => sb%next                     !point to next element in the list
   END DO
 END IF

 e => elset%head                       !point to first element
 iel = 0
 DO i=1,elset%nelem
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
   WRITE (50,ERR=9999) e%area1,e%lb,e%angle,e%a,e%b,e%c,e%cd,e%a0,e%ci,e%gamma, &
                       e%stra1,elset%stint(:,iel),e%mems
   IF( size1 > 0 )WRITE (50,ERR=9999) e%gausv
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)
 IF (elset%shear < 0 )THEN
   size1 = SIZE(elset%factors)
   WRITE(50,ERR=9999) size1
   WRITE(50,ERR=9999) (elset%factors(i),i=1,size1)
   WRITE(50,ERR=9999) (elset%ninv(i),i=1,npoin)
 END IF
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump13
