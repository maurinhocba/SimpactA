 SUBROUTINE dump15( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE
 TYPE (ele15_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs

 INTEGER (kind=4) :: i,size1,size2,iel
 TYPE (ele15), POINTER :: e
 TYPE (reinf), POINTER :: sb
 CHARACTER(len=midn):: straux

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'CBST'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, elset%nrf,  elset%logst, &
                     elset%lside, elset%gauss, elset%plstr, elset%angdf

 ! dump branching sides data base
 IF( elset%nrf > 0 )THEN
   sb => elset%rhead
   DO i=1,elset%nrf
     ! dump variables associated to the side
     WRITE (50,ERR=9999) sb%numel,sb%lnods,sb%secno
     WRITE (50,ERR=9999) sb%l0,sb%curv,sb%fc,sb%dc,sb%stint
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
                       e%stra1,elset%stint(:,iel) ,e%mems
   IF( size1 > 0 )WRITE (50,ERR=9999) e%gausv
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump15
