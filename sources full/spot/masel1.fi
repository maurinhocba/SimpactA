 SUBROUTINE masel1(ndime,nreqs,nnode,nelem,narch,heade,ngrqs,sname)
 IMPLICIT NONE
 INTEGER (kind=4) :: ndime,nreqs,nnode,nelem,narch
 INTEGER (kind=4), POINTER :: ngrqs(:)
 CHARACTER (len=*) :: sname
 TYPE (ele01), POINTER :: heade

 INTEGER (kind=4) ielem
 TYPE (ele01), POINTER :: e

 !     initial DATA for post-process

 WRITE(17,ERR=9) 1,nelem,sname   !1=SPOT
 !IF( nnode == 2 )THEN
   WRITE(17,ERR=9) nnode,3
 !ELSE
 !  WRITE(17,ERR=9) nnode,4
 !END IF
 e => heade
 DO ielem=1,nelem
   WRITE(17,ERR=9)e%matno,e%lnods
   IF( nnode == 6 ) WRITE(17,ERR=9)e%gausv(1:4)
   e => e%next
 END DO
 ! print auxiliar coordinates if necessary
 IF( nnode == 2 )THEN
   e => heade
   DO ielem=1,nelem
     IF( e%lnods(2) == 0 )WRITE(17,ERR=9)e%gausv(2:1+ndime)
     e => e%next
   END DO
 END IF
 !*** requested gauss points for output

 IF(nreqs > 0) THEN
   CALL openfi(nfile=21,narch=narch)
   IF( nnode == 2 )THEN
     ! 3 is the number of variables at gauss points stored for post-process
     WRITE(narch,ERR=9) nreqs,3,ngrqs(1:nreqs)            ! force,momen,epdef
     WRITE(narch,ERR=9) 1,sname,'N    ','M    ','EPD  '
   ELSE
     WRITE(narch,ERR=9) nreqs,6,ngrqs(1:nreqs)
     WRITE(narch,ERR=9) 1,sname,'N    ','T_1  ','T_2  ', &
                                   'R_1  ','R_2  ','R_3  '
   END IF
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE masel1
