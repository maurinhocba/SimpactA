 SUBROUTINE masel8(nnode,ngaus,nreqs,nelem,narch, &
                   heade,ngrqs,sname)
 IMPLICIT NONE
 INTEGER (kind=4), PARAMETER :: nstre = 6
 INTEGER (kind=4) :: nnode,ngaus,nreqs,nelem,narch
 INTEGER (kind=4), POINTER :: ngrqs(:)
 CHARACTER (len=*) :: sname
 TYPE (ele08), POINTER :: heade

 INTEGER (kind=4) n,i
 TYPE (ele08), POINTER :: e


 WRITE(17,ERR=9999) 8,nelem,sname
 WRITE(17,ERR=9999) nnode,nstre+2,ngaus
 e => heade
 DO i=1,nelem
   WRITE(17,ERR=9999) e%matno,(e%lnods(n),n=1,nnode)
   e => e%next
 END DO
 IF(nreqs > 0) THEN
   CALL openfi(nfile=21,narch=narch)
   WRITE(narch,ERR=9999) nreqs,nstre+2,(ngrqs(i),i=1,nreqs)
   WRITE(narch,ERR=9999) 8,sname,'N    ','Q_2  ','Q_3  ', &
        'M_T  ','M_2  ','M_3  ', 'M_T  ','EPD  ','Lambd'
 END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE masel8
