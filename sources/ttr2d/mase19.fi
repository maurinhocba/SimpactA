 SUBROUTINE mase19 (ntype,nreqs,nelem,heade,ngrqs,narch,angle,sname)

 ! stores mesh attributes for post-process

 IMPLICIT NONE

 INTEGER (kind=4) :: ntype,nreqs,nelem,narch
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle
 TYPE (ele19), POINTER  :: heade
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) i
 TYPE (ele19), POINTER :: e

 WRITE(17,ERR=9) 20,nelem,sname           !elm_type, number of elements, set_name
 WRITE(17,ERR=9) ntype,1,nnode,6          !ntype,ngaus,nnode,nstre
 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%matno,e%lnods(1:nnode)    !material, and nodes (main element)
   e => e%next                      !point to next element
 END DO

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%angle+angle           !local angle
   e => e%next                      !point to next element
 END DO

 IF(nreqs > 0) THEN                         !if output requested
   CALL openfi(nfile=21,narch=narch)        !open file
   WRITE(narch,ERR=9) nreqs,7,(ngrqs(i),i=1,nreqs) !write header
   SELECT CASE ( ntype )
   CASE (1,2)
     WRITE(narch,err=9) 19,sname,'S_xx ','S_yy ','S_xy ','S_zz ','Eps  ','Eqs  ','S_y  '
   CASE (3)
     WRITE(narch,err=9) 19,sname,'S_rr ','S_zz ','S_rz ','S_tt ','Eps  ','Eqs  ','S_y  '
   END SELECT
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE mase19
