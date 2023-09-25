 SUBROUTINE dump19( elset )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

   !--- Dummy variables
   TYPE(ele19_set),POINTER:: elset          !element set
   !--- Local variables
   INTEGER(kind=4)::  nreqs
   INTEGER (kind=4) :: i, nvare
   CHARACTER(len=midn):: straux
   TYPE (ele19), POINTER :: e

   straux = 'ELSET'
   WRITE(50,ERR=9999) straux               !key word
   straux = 'TR2D'
   WRITE(50,ERR=9999) straux,elset%sname    !element type and element set name
   WRITE(50,ERR=9999) elset%nelem           !number of elements
   WRITE(50,ERR=9999) elset%nelem, elset%nreqs, elset%narch, &
                      elset%gauss, elset%lside, elset%eulrf, elset%plstr, elset%angdf
   WRITE(50,ERR=9999) (elset%eside(:,i),i=1,elset%nelem)
   ! dump element numbers where output is required
   nreqs = elset%nreqs
   IF (nreqs > 0) WRITE(50,ERR=9999) (elset%ngrqs(i), i=1,nreqs)

   e => elset%head                       !point to first element
   DO
     IF (.NOT.ASSOCIATED(e)) EXIT  !last element processed, EXIT loop
     ! dump variables associated to an element
     IF( ASSOCIATED(e%gausv) )THEN
       nvare = SIZE (e%gausv)
     ELSE
       nvare = 0
     END IF
     WRITE(50,ERR=9999) e%numel,e%matno,e%lnods, nvare
     WRITE(50,ERR=9999) e%area1,e%angle,e%cd,e%stint,e%beta,e%repla
     IF( nvare > 0 )WRITE (50,ERR=9999) e%gausv
     e => e%next                     !point to next element in the list
   END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump19
