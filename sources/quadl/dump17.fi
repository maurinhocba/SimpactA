 SUBROUTINE dump17( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

 TYPE (ele17_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs

 INTEGER (kind=4) :: i,j,nucom,nvare, osec
 CHARACTER(len=midn):: straux
 TYPE (ele17), POINTER :: e
 TYPE (section), POINTER :: sec
 LOGICAL :: found,plast

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'QUADL'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, &
            elset%gauss, elset%plstr, &
            elset%angdf, elset%ngaus
 osec = -1
 e => elset%head                       !point to first element
 nvare = 5           !number of internal variables per component
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT  !last element processed, EXIT loop
   IF( e%matno /= osec )THEN
     sec => psecs(e%matno)%p 
     nucom = sec%iprop(2)       !number of components
     plast = sec%mtbas%matdef(3) > 1
     osec = e%matno
   END IF
   ! dump variables associated to an element
   WRITE (50,ERR=9999) e%numel,e%matno,e%lnods
   WRITE (50,ERR=9999) e%dvol,e%angle,e%cartd,e%stint
   IF( plast ) &
     WRITE(50,ERR=9999) ((e%gausv(i,j),i=1,nvare*nucom),j=1,elset%ngaus) !Write internal variables
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump17
