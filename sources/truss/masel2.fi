 SUBROUTINE masel2(nreqs,nelem,narch,heade,ngrqs,sname)
 IMPLICIT NONE
 INTEGER (kind=4) :: nreqs,nelem,narch
 INTEGER (kind=4), POINTER :: ngrqs(:)
 CHARACTER (len=*) :: sname
 TYPE (ele02), POINTER :: heade

 INTEGER (kind=4) ielem
 TYPE (ele02), POINTER :: e

 !     initial DATA for post-process

 WRITE(17,ERR=9) 2,nelem,sname
 WRITE(17,ERR=9) nnode,3
 e => heade
 DO ielem=1,nelem
   WRITE(17,ERR=9) e%matno,e%lnods
   e => e%next
 END DO

 !*** requested gauss points for output

 IF(nreqs > 0) THEN
   CALL openfi(nfile=21,narch=narch)
   WRITE(narch,ERR=9) nreqs,3,ngrqs(1:nreqs)
   WRITE(narch,ERR=9) 2,sname,'N    ','Str  ','EPD  '
 END IF
 !     3 is the number of variables at gauss points stored for post-process

 RETURN
 9 CALL runen2('')
 END SUBROUTINE masel2
