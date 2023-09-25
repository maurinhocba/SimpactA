 SUBROUTINE inpco3
!***********************************************************************
!
! **  Modify contact interfaces parameters at a new Strategy
!     slave and master surfaces can not be changed
!
!***********************************************************************
 IMPLICIT NONE

!     Local Variables
 INTEGER (kind=4) i
 TYPE (pair3_db), POINTER :: pair
 CHARACTER(len=mnam):: pname
 LOGICAL :: found

 WRITE(lures,"(//,'     Modified Contact Parameters',/)",ERR=9999)
 CALL listen('INPCO3')                        !read PAIR_INFORMATION card
 IF(.NOT.exists('PAIRIN'))CALL runend &
   ('INPCO3: PAIR_INFO card expected !!!')
 DO
   CALL listen('INPCO3')                      !read pair parameters
   IF( exists('ENDPAI')) EXIT                 !exit if END_PAIR card read

   pname = get_name('IPAIR ',found,'! PAIR to be modified ',stype='PAIR') !Pair name

   ! search pair to modify
   i = 0                                      !initializes
   pair => headp                              !point to first pair
   DO
     i=i+1
     IF(TRIM(pair%pname) == TRIM(pname))THEN
       EXIT
     ELSE IF(i == npair)THEN
       WRITE(lures,"(' pair ',A,' DOES NOT exist ')",ERR=9999) TRIM(pname)
       CALL runend('INPCO3: Non-existing pair to modify')
     END IF
     pair => pair%next
   END DO
!....   read NEW contact surfaces properties data for the pair
   pair%indcon=getint('INDCON',pair%indcon, &
                      ' Contact formulation flag .........')
   WRITE(lures,"(12x,'0 flexible-flexible ',/12x,'1 flexible-rigid    ')",ERR=9999)
   pair%npenal=getrea('NPENAL',pair%npenal, &
                      ' Normal penalty coeff .............')
   pair%tpenal=getrea('TPENAL',pair%tpenal, &
                      ' Tangential penalty coeff..........')
   pair%static=getrea('STATIC',pair%static, &
                      ' Static Friction coefficient ......')
   IF( exists('STATIC') )pair%kinet = pair%static
   pair%kinet =getrea('KINET ',pair%kinet , &
                      ' Kinetic Friction coefficient .....')
   IF( therm )THEN
     pair%hcont =getrea('HCONT ',pair%hcont, &
                ' H factor for heat conduction......')
     pair%econt =getrea('EHCON ',pair%econt, &
                ' Exponent for heat conduction......')
     pair%hv    =getrea('VICKH ',pair%hv,    &
                ' Vickers  hardness ................')
     pair%rec   =getrea('EFFUS ',pair%rec,   &
                ' Relative Effusivity (Slave).......')
   END IF
   pair%start =getrea('START ',pair%start , &
                      ' Activation time ..................')
   pair%end   =getrea('END   ',pair%end   , &
                      ' Deactivation time  ...............')
   pair%cutoff=-getrea('CUTOFF',-pair%cutoff,   &
                      ' CUT OFF value for penetration.....')
   IF( pair%cutoff >= 0d0 )pair%cutoff = -1d9
   pair%gapinc=-getrea('GAPINC',-pair%gapinc,      &
                      ' Maximum Incremental Gap per Step..')
   IF( pair%gapinc >= 0d0 ) pair%gapinc = pair%cutoff/20
   pair%freq  = getint('FREQ  ',pair%freq,  &
                  ' Frequency of global search .......')
   IF( exists('PRESS ') )THEN
     pair%press = .TRUE.
   ELSE IF(pair%press .AND. exists('NOPRES') )THEN
     pair%press = .FALSE.
     DEALLOCATE( pair%presn )
   END IF
 END DO

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE inpco3
