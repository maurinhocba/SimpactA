 SUBROUTINE mase13 (nreqs,nelem,heade,ngrqs,narch,angle,sname,locax)

 ! Prints initial data for postprocess

 IMPLICIT NONE
 INTEGER (kind=4) :: nreqs,nelem,narch,locax
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle
 TYPE (ele13), POINTER  :: heade
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) i,lnods(6)
 TYPE (ele13), POINTER :: e
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (postv), POINTER :: postp       !pointer to first postprocess variable

 WRITE(17,ERR=9) 14,nelem,sname     !elm_type, number of elements, set_name
 e => heade                         !point to first element
 sec => psecs(e%matno)%p            !point to section of first element in the set
 IF( sec%secdef(4) > 0 )THEN        !if special post-process variables
   WRITE(17,ERR=9) 3,1,13+sec%secdef(5),locax     !nnode,ngaus,nstre,locax
   postp => sec%postp
   WRITE (17,ERR=9) sec%secdef(4),e%matno  !writes position of section in list
   DO i=1,sec%secdef(4)
     WRITE (17,ERR=9) postp%type,postp%dim,postp%name
     postp => postp%next
   END DO
 ELSE
   WRITE(17,ERR=9) 3,1,13,locax             !nnode,ngaus,nstre,locax
 END IF
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   lnods = e%lnods
   DO i=1,6
     IF(lnods(i) < 0 )lnods(i) = 0
   END DO
   WRITE(17,ERR=9) e%matno,lnods(1:6)    !material, and nodes (whole patch)
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
   WRITE(narch,ERR=9) nreqs,8,(ngrqs(i),i=1,nreqs) !write header
   WRITE(narch,ERR=9) 13,sname,'N_xx ','N_yy ','N_xy ','M_xx ','M_yy ','M_xy ', &
                'Q_x  ','Q_y  '
 END IF

 RETURN
 9 CALL runen2(' error while writing to the disk')
 END SUBROUTINE mase13
