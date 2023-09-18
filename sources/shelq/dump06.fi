 SUBROUTINE dump06( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE
 TYPE (ele06_set), POINTER :: elset          !element set

 INTEGER (kind=4) ::  nreqs

 INTEGER (kind=4) :: i
 CHARACTER(len=midn):: straux
 TYPE (ele06), POINTER :: e

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux                !key word
 straux = 'SHELQ'
 WRITE (50,ERR=9999) straux,elset%sname    !element type and element set name
 WRITE (50,ERR=9999) elset%nelem           !number of elements
 WRITE (50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, elset%gauss, &
                     elset%plstr, elset%angdf, elset%zigzag, elset%nstre

 e => elset%head                       !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT  !last element processed, EXIT loop
   ! dump variables associated to an element
   WRITE (50,ERR=9999) e%numel,e%matno,e%lnods
   WRITE (50,ERR=9999) e%angle,e%dvolu,e%cartd,e%jacin,e%stra0,e%gamm0,e%strsg,e%ambda
   IF( ASSOCIATED( e%ehist ) ) WRITE(50,ERR=9999) e%ehist,e%strap
   IF( ASSOCIATED( e%stres ) ) WRITE(50,ERR=9999) e%stres
   e => e%next                     !point to next element in the list
 END DO
 ! dump element numbers where output is required
 nreqs = elset%nreqs
 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump06
