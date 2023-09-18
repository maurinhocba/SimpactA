 SUBROUTINE arisi0 (naris,iwrit,actio)

 ! reads data for nodes fixed on a side (narisn)

 USE c_input
 USE nar_db
 IMPLICIT NONE
 CHARACTER(len=*),INTENT(INOUT):: actio
 INTEGER(kind=4):: naris,iwrit

 ! Local
 TYPE (nar_nod), POINTER :: nar,head,tail
 INTEGER :: n
 LOGICAL :: slide

 IF (iwrit == 1) WRITE(lures,"(//,3x,'Slave Nodes Data',/)",ERR=9999)

 !Initialize empty list
 CALL ini_nar(head,tail)
 IF (.NOT.exists('ADD   ')) THEN
   naris=0
 ELSE IF (exists('ADD   ')) THEN

   IF (iwrit == 1) WRITE(lures,"(/,3x, &
&   'Kept from the previous strategy nodes slaved to arista.',/)",ERR=9999)

   DO n=1,naris                 !transfer data from array to list
     ALLOCATE (nar)             !get memory
     nar%slave  = nardat(1,n)   !transfer present data
     nar%mast01 = nardat(2,n)
     nar%mast02 = nardat(3,n)
     CALL add_nar (nar, head, tail)  !add to list
   END DO
 END IF

 IF (TRIM(actio) == 'NSTRA0') actio='NSTRA1'
 IF (ASSOCIATED (nardat) ) DEALLOCATE (nardat)  !release memory

 IF (iwrit == 1) WRITE(lures,"(/'   Slave Node       Side Nodes   Slide')",ERR=9999)

 !Loop to read data and add them to the list
 DO
   CALL listen('ARISI0')               !read a card
   IF (exists('ENDNAR')) EXIT          !if end of data found, exit
   ALLOCATE (nar)                      !get memory
   IF( exists('SLAVE') )THEN
     nar%slave  = getint('SLAVE ',0,'!Number of a Slave node ...........') ! slave node
     nar%mast01 = getint('MAST01',0,'!Number of first Master node ......') ! first master node
     nar%mast02 = getint('MAST02',0,'!Number of second  Master node ....') ! second master node
     slide = exists('SLIDE')
   ELSE
     nar%slave  = INT(param(1)) ! slave node
     nar%mast01 = INT(param(2)) ! first master node
     nar%mast02 = INT(param(3)) ! second master node
     slide      = INT(param(4)) == 1
   END IF
   naris=naris+1                       !increase counter
   IF (iwrit == 1) WRITE(lures, '(3i10,l8)',ERR=9999)         & !echo data (twice)
           nar%slave, nar%mast01, nar%mast02, slide
   IF( slide )nar%slave = -nar%slave   !sliding node
   CALL add_nar( nar, head, tail )     !add to list
 END DO

 !Store data in an array
 ALLOCATE (nardat(3,naris))    !reserve space
 CALL store_nar(head)          !transfer data and release list space
 CALL ini_nar(head,tail)       !nullify pointers

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE arisi0
