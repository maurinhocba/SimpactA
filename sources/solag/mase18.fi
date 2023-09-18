 SUBROUTINE mase18 (nreqs,nelem,heade,ngrqs,narch,angle,sname,ngaus,locax)

 !performs the final tasks in ACVD18 at a restart

 IMPLICIT NONE

 INTEGER (kind=4) :: nreqs,nelem,narch,ngaus,locax
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle(:)
 TYPE (ele18), POINTER  :: heade
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) i
 TYPE (ele18), POINTER :: e
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (postv), POINTER :: postp       !pointer to first postprocess variable

 e => heade                         !point to first element
 WRITE(17,ERR=9) 18,nelem,sname  !elm_type, number of elements, set_name
 sec => psecs(e%matno)%p            !point to section of first element in the set
 IF( sec%secdef(4) > 0 )THEN        !if special post-process variables
   WRITE(17,ERR=9) ngaus,nnode,8+sec%secdef(5),locax
   postp => sec%postp
   WRITE (17,ERR=9) sec%secdef(4),e%matno  !writes position of section in list
   DO i=1,sec%secdef(4)
     WRITE (17,ERR=9) postp%type,postp%dim,postp%name
     postp => postp%next
   END DO
 ELSE
   WRITE(17,ERR=9) ngaus,nnode,8,locax  ! 8=nvarg
 END IF
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%matno,e%lnods(:)  !material, and nodes
   e => e%next                      !point to next element
 END DO

 WRITE(17,ERR=9) angle         !set angles
 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%angle       !local angles
   e => e%next                      !point to next element
 END DO

 IF(nreqs > 0) THEN                         !if output requested
   CALL openfi(nfile=21,narch=narch)        !open file
   WRITE(narch,ERR=9)nreqs,9,(ngrqs(i),i=1,nreqs) !write header
   WRITE(narch,ERR=9) 18,sname,'S_xx ','S_yy ','S_zz ','S_xy ','S_xz ', &
                               'S_yz ','Eps  ','Eqs  ','S_y  '
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE mase18
