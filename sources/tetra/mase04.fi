 SUBROUTINE mase04 (nreqs,nelem,heade,ngrqs,narch,angle,sname)

 !performs the final tasks in ACVD04 at a restart

 IMPLICIT NONE
 INTEGER (kind=4) :: nreqs,nelem,narch
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle(:)
 TYPE (ele04), POINTER  :: heade
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) i
 TYPE (ele04), POINTER :: e
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (postv), POINTER :: postp       !pointer to first postprocess variable

 e => heade                         !point to first element
 WRITE(17,ERR=9) 5,nelem,sname      !elm_type, number of elements, set_name
 sec => psecs(e%matno)%p            !point to section of first element in the set
 IF( sec%secdef(4) > 0 )THEN        !if special post-process variables
   WRITE(17,err=9) 1,4,8+sec%secdef(5),0 !ngaus,nnode,nstre,locax
   postp => sec%postp
   WRITE (17,err=9) sec%secdef(4)
   DO i=1,sec%secdef(4)
     WRITE (17,err=9) postp%type,postp%dim,postp%name
     postp => postp%next
   END DO
 ELSE
   WRITE(17,ERR=9) 1,4,8,0            !1=ngaus for each direction 4=nnode 8=nvarg
 END IF
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%matno,e%lnods(1:4)  !material, and nodes
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
   WRITE(narch,err=9) 4,sname,'S_xx ','S_yy ','S_zz ','S_xy ','S_xz ', &
                              'S_yz ','Eps  ','Eqs  ','S_y  '
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE mase04
