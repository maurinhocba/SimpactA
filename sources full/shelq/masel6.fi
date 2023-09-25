 SUBROUTINE masel6 (nreqs,nelem,heade,ngrqs,narch,angle,sname,nstre)

 !performs the final tasks in ACVDF6 at a restart

 IMPLICIT NONE
 INTEGER (kind=4) :: nreqs,nelem,narch,nstre
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle
 TYPE (ele06), POINTER  :: heade
 CHARACTER (len=*) :: sname
 !
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (postv), POINTER :: postp       !pointer to first postprocess variable
 INTEGER (kind=4) i
 TYPE (ele06), POINTER :: e

 WRITE(17,ERR=9) 6,nelem,sname           !elm_type, number of elements, set_name

 e => heade                        !point to first element
 sec => psecs(e%matno)%p            !point to section of first element in the set
 IF( nstre == 14 .AND. sec%secdef(4) > 0 )THEN        !if special post-process variables
   WRITE(17,ERR=9) nnode,4,11+sec%secdef(5)           !nnode,ngaus,nstre
   postp => sec%postp
   WRITE (17,ERR=9) sec%secdef(4),e%matno
   DO i=1,sec%secdef(4)
     WRITE (17,ERR=9) postp%type,postp%dim,postp%name
     postp => postp%next
   END DO
 ELSE
   WRITE(17,ERR=9) nnode,ngaus,11  !nnode,ngaus,nstre+3
 END IF

 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%matno,e%lnods(1:nnode)  !material, and nodes
   e => e%next                      !point to next element
 END DO

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,ERR=9) e%angle+angle  !local angle
   e => e%next                      !point to next element
 END DO

 IF(nreqs > 0) THEN                         !if output requested
   CALL openfi(nfile=21,narch=narch)        !open file
   WRITE(narch,ERR=9) nreqs,11,(ngrqs(i),i=1,nreqs) !write header
   WRITE(narch,err=9) 6,sname,'N_xx ','N_yy ','N_xy ','M_xx ','M_yy ','M_xy ', &
                'Q_x  ','Q_y  ','EPS_B','EPS_T','Th_rt'
 END IF

 RETURN
 9 CALL runen2('')
 END SUBROUTINE masel6
