 SUBROUTINE inpco1()
!***********************************************************************
! Modify contact interfaces parameters at a new Strategy
! slave and master surfaces can not be changed
!***********************************************************************
 IMPLICIT NONE

 !--- Local Variables
 INTEGER(kind=4):: i
 TYPE(pair1_db),POINTER:: pair
 CHARACTER(len=mnam):: pname
 LOGICAL :: found

 WRITE(lures,"(//,'     Modified Contact Parameters',/)",ERR=9999)
 CALL listen('INPCO1')                        !read PAIR_INFORMATION card
 IF(.NOT.exists('PAIRIN'))CALL runend('INPCO1: PAIR_INFO card expected !!!')
 DO
   CALL listen('INPCO1')                      !read pair parameters
   IF (exists('ENDPAI')) EXIT                 !exit if END_PAIR card read

   pname = get_name('IPAIR ',found,'! PAIR to be modified ',stype='PAIR') !Pair name

   ! search pair to modify
   i = 0                                      !initializes
   pair => headp                              !point to first pair
   DO
     i = i + 1
     IF (TRIM(pair%pname) == TRIM(pname)) THEN
       EXIT
     ELSE IF(i == npair)THEN
       WRITE(lures,"(' pair ',A,' DOES NOT exist ')",ERR=9999) TRIM(pname)
       CALL runend('INPCO1: Non-existing pair to modify')
     END IF
     pair => pair%next
   END DO
   !....   read NEW contact surfaces properties data for the pair
   pair%indcon = getint('INDCON',pair%indcon,' Contact formulation flag .........')
   WRITE(lures,"(12x,'0 flexible-flexible ',/12x,'1 flexible-rigid    ')",ERR=9999)
   pair%npenal = getrea('NPENAL',pair%npenal,' Normal penalty coeff .............')
   pair%tpenal = getrea('TPENAL',pair%tpenal,' Tangential penalty coeff..........')
   pair%static = getrea('STATIC',pair%static,' Static Friction coefficient ......')
   IF (exists('STATIC')) pair%kinet = pair%static
   pair%kinet  = getrea('KINET ',pair%static,' Kinetic Friction coefficient .....')
   pair%start  = getrea('START ',pair%start ,' Activation time ..................')
   pair%end    = getrea('END   ',pair%end   ,' Deactivation time  ...............')
   IF (exists('PRESS ')) THEN
     pair%press = .TRUE.
   ELSE IF (pair%press .AND. exists('NOPRES')) THEN
     pair%press = .FALSE.
     DEALLOCATE(pair%presn)
   END IF
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpco1
