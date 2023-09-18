 SUBROUTINE mase17 (ntype,nreqs,nelem,heade,ngrqs,narch,angle,sname,ngaus)

 !performs the final tasks in ACVD17 at a restart

 IMPLICIT NONE

 INTEGER (kind=4) :: ntype,nreqs,nelem,narch,ngaus
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle
 TYPE (ele17), POINTER  :: heade
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) i,ngaud!,isec,lblflc
 TYPE (ele17), POINTER :: e
 ! TYPE(section),POINTER:: sec    !pointer to a section data

 !isec = heade%matno                !associated section
 !sec => psecs(isec)%p              !point to section of first element in the set
 !lblflc = sec%iprop(3)             !FLC curve label
 !ngaud = NINT(SQRT(REAL(ngaus)))
 ngaud = 2
 IF(ngaus == 1)ngaud=1
 WRITE(17,ERR=9) 17,nelem,sname           !elm_type, number of elements, set_name
 WRITE(17,ERR=9) ntype,ngaud,nnode,6      !6=nstre
 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT          !if all elements processed, exit loop
   WRITE(17,ERR=9) e%matno,e%lnods(:)   !material, and nodes
   e => e%next                             !point to next element
 END DO

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%angle+angle !local angle
   e => e%next                      !point to next element
 END DO

 IF(nreqs > 0) THEN                         !if output requested
   CALL openfi(nfile=21,narch=narch)        !open file
   WRITE(narch,ERR=9) nreqs,7,(ngrqs(i),i=1,nreqs) !write header
   SELECT CASE ( ntype )
   CASE (1,2)
     WRITE(narch,err=9) 17,sname,'S_xx ','S_yy ','S_xy ','S_zz ','Eps  ','Eqs  ','S_y  '
   CASE (3)
     WRITE(narch,err=9) 17,sname,'S_rr ','S_zz ','S_rz ','S_tt ','Eps  ','Eqs  ','S_y  '
   END SELECT
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE mase17
