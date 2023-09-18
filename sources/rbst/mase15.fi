 SUBROUTINE mase15 (nreqs,nelem,heade,ngrqs,narch,angle,sname,nrf,rhead)

 ! Prints initial data for postprocess

 IMPLICIT NONE
 INTEGER (kind=4) :: nreqs,nelem,narch,nrf
 INTEGER (kind=4), POINTER  :: ngrqs(:)
 REAL (kind=8 ) :: angle
 TYPE (ele15), POINTER  :: heade
 TYPE (reinf), POINTER  :: rhead
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) i,lnods(6)
 TYPE (ele15), POINTER :: e
 TYPE (reinf), POINTER :: si
 CHARACTER (len=mnam) :: rname
 TYPE (section), POINTER :: sec    !pointer to a section data
 TYPE (postv), POINTER :: postp       !pointer to first postprocess variable

 WRITE(17,err=9) 14,nelem,sname     !elm_type, number of elements, set_name
 e => heade                         !point to first element
 sec => psecs(e%matno)%p            !point to section of first element in the set
 IF( sec%secdef(4) > 0 )THEN        !if special post-process variables
   WRITE(17,err=9) 3,1,13+sec%secdef(5),3   !nnode,ngaus,nstre,locax
   postp => sec%postp
   WRITE (17,err=9) sec%secdef(4)
   DO i=1,sec%secdef(4)
     WRITE (17,err=9) postp%type,postp%dim,postp%name
     postp => postp%next
   END DO
 ELSE
   WRITE(17,err=9) 3,1,13,3                 !nnode,ngaus,nstre,locax
 END IF
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   lnods = e%lnods
   DO i=1,6
     IF(lnods(i) < 0 )lnods(i) = 0
   END DO
   WRITE(17,err=9) e%matno,lnods(1:6)    !material, and nodes (whole patch)
   e => e%next                             !point to next element
 END DO

 e => heade                         !point to first element
 DO
   IF (.NOT.ASSOCIATED (e) ) EXIT   !if all elements processed, exit loop
   WRITE(17,err=9) e%angle+angle           !local angle
   e => e%next                      !point to next element
 END DO

 IF(nreqs > 0) THEN                         !if output requested
   CALL openfi(nfile=21,narch=narch)        !open file
   WRITE(narch,err=9) nreqs,8,(ngrqs(i),i=1,nreqs) !write header
   WRITE(narch,err=9) 15,sname,'N_xx ','N_yy ','N_xy ','M_xx ','M_yy ','M_xy ', &
                'Q_x  ','Q_y  '
 END IF

 IF( nrf > 0 )THEN
   rname = sname
   i = LEN_TRIM(sname)
   rname(i+1:i+3) = '_rf'
   WRITE(17,err=9)  9,nrf,rname
   si => rhead                         !point to first element
   WRITE(17,err=9) 2,6,1,5             !nnode,ngaus,nstre,ntype = 5 is a flag
   DO
     IF (.NOT.ASSOCIATED (si) ) EXIT   !if all elements processed, exit loop
     WRITE(17,err=9) si%secno,si%lnods(1:2)    !material, and nodes (whole patch)
     si => si%next                      !point to next element
   END DO

 END IF

 RETURN
 9 CALL runen2(' error while writing to the disk')
 END SUBROUTINE mase15
