 SUBROUTINE dumpi8(sname, nelem, nnode, ngaus, axesc,         &
                   nreqs, narch, heade,                &
                   ngrqs, posgp, shape, deriv, weigh)
 !*****************************************************************
 !
 !*** dumps an element set for restart
 !
 !*****************************************************************
 IMPLICIT NONE

 CHARACTER(len=*):: sname
 INTEGER(kind=4):: nelem, nnode, ngaus, axesc, nreqs, narch
 INTEGER (kind=4), POINTER :: ngrqs(:)
 REAL (kind=8) :: posgp(:), shape(:,:), deriv(:,:), weigh(:)
 TYPE (ele08), POINTER :: heade

 INTEGER (kind=4) :: i,ielem
 CHARACTER(len=midn):: straux
 TYPE (ele08), POINTER :: e

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux
 straux = 'BEAM'
 WRITE (50,ERR=9999) straux,sname
 WRITE (50,ERR=9999) nelem
 WRITE (50,ERR=9999) nelem, nnode, ngaus, axesc, nreqs, narch

 WRITE (50,ERR=9999) posgp, shape, deriv, weigh
 e => heade
 DO ielem=1,nelem
   WRITE (50,ERR=9999) e%lnods(1:nnode) , e%matno, e%numel
   WRITE (50,ERR=9999) e%llbd,e%jac,e%stra0,e%stran,e%stres,e%epdef,e%sedef
   e => e%next
 END DO

 IF ( nreqs > 0 ) WRITE (50,ERR=9999) (ngrqs(i), i=1,nreqs)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dumpi8
