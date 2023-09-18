 SUBROUTINE beginp(actio,nrest,cpui,kstep)
!*********************************************************************
!
!     reads configuration file
!           state if its a NEW problem or a RESTART, etc
!           input data file name
!           output file name-root
!           restart file name,
!           time or steps between dumping of restar file
!           steps between dumping screen updates
!           memory size for auxiliar array
!
!*********************************************************************
 USE param_db,ONLY: mlen, mprg, max_memo
 USE name_db,ONLY : input, output, namr, prognm, data_file, gen_data_file
 USE c_input,ONLY : openfi, del_old_files, upcase
 USE kinc_db, ONLY : nn
 USE ctrl_db, ONLY : memo
 USE gvar_db, ONLY : actchk

 IMPLICIT NONE

   ! dummy arguments
   CHARACTER(len=*),INTENT (OUT) :: actio     ! NEW, RESTAR or NSTRA0
   INTEGER (kind=4),INTENT (OUT) :: nrest,  & ! time between dumping of restar files [sec]
                                    kstep     ! steps between screen updating
   REAL(kind=8),INTENT(OUT) :: cpui           ! CPU used for present model

   ! local variables
   CHARACTER(len=1),PARAMETER :: cn = '@'  ! character no to generate data_file
   INTEGER(kind=4),PARAMETER :: npar = 9,    & !Different input number of program parameters
                                mpar = 4       !Maximum length for input commands
   INTEGER(kind=4) :: narg,i,ichek,lengs(npar),numpar(npar),j,k,number
   LOGICAL :: error(npar),errora
   CHARACTER(len=mpar),PARAMETER :: words(npar)=(/'    ','INP=','OUT=','FRQ=','RSF=','MEM=',   &
                                                  'FRS=','LST ','CHK '/)
   CHARACTER(len=2*mlen) :: param0       !auxiliar string
   CHARACTER(len=mlen) :: params(npar)   !parameters
   CHARACTER (len=mlen) :: straux        !auxiliar string
   CHARACTER(len=mpar) :: word           !auxiliar short string


   INTERFACE
     INCLUDE 'wrtrpt.h'
   END INTERFACE
   INTEGER(kind=4) :: nargs,  &  !number of arguments in command line (system routine)
                      istat      !input error flag

   narg = nargs()                !number of arguments in command line
   narg = narg-1                 !number of optional parameters (exclude program itself)

!-------------------set defaults
   actio = 'NEW'    !new model
   input = ''       !input data file name
   output= ''       !root name for output data files
   nrest = 1800     !half an hour in seconds
   namr  = ''       !name of the re-start file
   memo  = 13*nn    !size (integer) of the auxiliar array
   kstep = 1250     !number of steps between screen updates
   error = .TRUE.   !initializae flag array
! for optional arguments initializes flags
   error(1) = .FALSE.  ! actio
   error(4) = .FALSE.  ! frq
   error(6:9) = .FALSE.  ! mem, frs, lst, chk

   param0 = ''           !initializes auxiliar string
!------------------- IF no arguments use COMAN.DAT file
   IF (narg < 1) THEN      !no optional arguments in command line
     OPEN(1,FILE='coman.dat',FORM='formatted',STATUS='old',IOSTAT=ichek)
     IF (ichek /= 0) THEN  !coman.dat file does not exist in present directory
       WRITE(*,1000)       !write help information
       STOP                !and stop
     END IF
     READ(1,'(A)',IOSTAT=ichek) param0   !read first line
     i=1                   !initializes counter
     DO                    !loop until end of the file or i > npar
       lengs(i) = LEN_TRIM(param0)  !length of line read
       IF( lengs(i) > 0 )THEN
         params(i) = param0
         i = i+1             !update counter
       END IF
       IF (i > npar) EXIT   !EXIT if maximum number of params read
       READ(1,'(a)',IOSTAT=ichek) param0   !read a line
       IF (ichek /= 0) EXIT   !if read error exit
     END DO
     narg = i-1            !update number of arguments read successfully
     CLOSE(1)              !close Configuration file

   ELSE IF(narg > 1) THEN  !-------------------extract all arguments from command line
     params(1:npar) = ''  !initializes
     DO i=1,narg            !for each optional argument
       CALL getarg(i,params(i),istat)                      !intel_fc
       lengs(i) = LEN_TRIM(params(i)) !length of parameter read
     END DO
   ELSE                    !Only one argument, it must be the Input data file
     CALL getarg(1,param0,istat)                      !PC
     params(1) = param0
     lengs(1) = LEN_TRIM(params(1)) !length of parameter read
   END IF


   IF (narg == 1) THEN
   !-------------------assume datafile name IF only one argument
     word(1:mpar) = params(1)(1:mpar)     !first 4 chars
     CALL upcase(word)                    !upper case version
     !  assign name to INPUT
     IF (word(1:mpar) == 'INP=') THEN     !check if 'INP=' was included
       input(1:lengs(1)-mpar) = params(1)(mpar+1:lengs(1))  !remove 'INP='
     ELSE
       input(1:lengs(1)) = params(1)(1:lengs(1))
     END IF
     !            create output file name
     IF( input(1:1) == cn  )THEN
       gen_data_file = .FALSE.
       input = input(2:)
     END IF
     j = SCAN(input,'.', BACK = .TRUE.) -1
     IF( j > 0 )THEN
       output = input(1:j)
     ELSE
       output = input
     END IF
     error(2:3) = .FALSE.
     error(5)   = .FALSE.
     param0 = TRIM(input)
     CALL upcase(param0)
     IF (INDEX(param0,'.DAT')==0) THEN
       WRITE(*,1000)
       STOP
     END IF
   ELSE      !-------------------more than one argument
     DO i=1,narg         !for each argument read
       numpar(i) = 0     !initializes number of parameters
       word(1:mpar)=params(i)(1:mpar)  !read into WORD
       CALL upcase(word)               !upcase WORD
       ! check if main key-word is present
       !  first: check Main keyword
       IF (word(1:mpar)=='NEW ' .OR. word(1:mpar)=='RES ' .OR. word(1:mpar)=='NST ') THEN
         numpar(i) = 1           !first
         CALL upcase(params(i))
         CYCLE
       END IF
       ! second: check other keywords
       IF ((INDEX(params(i),'=') == 0) .OR. (lengs(i) <= mpar)) THEN !no parameter
         errora = .TRUE.   !         !  assume that key-word INP is not present
         IF (INDEX(params(i),'.') /= 0) THEN  !if a file name with extension is present
           errora = .FALSE.                   !no error
           numpar(i) = 100                    !flag to read input file
         END IF
         ! check if it is not a keyword without parameter
         IF (TRIM(word)==TRIM(words(8)).OR.TRIM(word)==TRIM(words(9))) errora=.FALSE.
         IF (errora) THEN  !invalid input, stop
           WRITE(*,1000)
           STOP
         END IF
       END IF
       DO j=2,npar  !check the existence of keywords
         ichek = INDEX(word,words(j))
         IF(ichek == 1) THEN
           numpar(i) = j
           EXIT
         END IF
       END DO
     END DO
     !--------------------
     DO i=1,narg

       SELECT CASE (numpar(i))
       CASE (1)               ! actio
                              ! if not input default 'NEW'
         IF(INDEX(params(i),'NEW') == 1) THEN ! "NEW"    : start a new problem
           actio = 'NEW'
         ELSE IF (INDEX(params(i),'RES') == 1) THEN  ! "RESTAR" : restart an old problem
           actio = 'RESTAR'                         !            without any modifications
           IF (error(8)) STOP 'ERROR: "RES" and "LST" can''t be defined simultaneously.'
         ELSE IF (INDEX(params(i),'NST') == 1) THEN  ! "NSTRA0" : restart an old problem
           actio = 'NSTRA0'                         !            and start a new strategy
         END IF
       CASE (2)               ! input file
         input(1:lengs(i)-mpar) = params(i)(mpar+1:lengs(i))
         IF( input(1:1) == cn )THEN
           gen_data_file = .FALSE.
           input = input(2:)
         END IF
         error(2) = .FALSE.

       CASE (3)               ! output file
         output = params(i)(mpar+1:lengs(i))
         error(3) = .FALSE.

       CASE (4)               ! restart dumping frequency
         k = lengs(i) - mpar
         straux(1:k) = params(i)(mpar+1:lengs(i))
         straux(k+1:k+1) = ' '
         nrest = number(straux)              !read in minutes
         IF (nrest > 0)   nrest = nrest*60   !store in seconds
         ! if nrest < 0    is the number of steps between dumping
       CASE (5)               ! restart file name
         namr(1:lengs(i)-mpar) = params(i)(mpar+1:lengs(i))
         error(5) = .FALSE.

       CASE (6)               ! memory for auxiliar array
         k = lengs(i) - mpar
         straux(1:k) = params(i)(mpar+1:lengs(i))
         straux(k+1:k+1) = ' '
         memo = number(straux)
         IF (memo < 256*1024) memo=memo*1024

       CASE (7)               ! number of steps between screen updates
         k = lengs(i) - mpar
         straux(1:k) = params(i)(mpar+1:lengs(i))
         straux(k+1:k+1) = ' '
         kstep = number(straux)

       CASE (8)               ! Restart from LaST restart file
         error(8) = .TRUE.
         IF (TRIM(actio)=='RESTAR') STOP 'ERROR: "RES" and "LST" can''t be defined simultaneously.'

       CASE (9)               ! check data file
         actchk = .TRUE.

       CASE (100)             ! input file without kwd
         input(1:lengs(i)) = params(i)(1:lengs(i))
         error(2) = .FALSE.

       END SELECT

     END DO
     !                    use input filename
     IF(error(3))THEN
       j = SCAN(input,'.', BACK = .TRUE.) -1
       IF( j > 0 )THEN
         output = input(1:j)
       ELSE
         output = input
       END IF
       error(3) = .FALSE.
     END IF

     IF (error(8)) THEN    !restart from last restart file
       namr = TRIM(output)
       CALL openfi(nfile=63,rest=namr)
       IF (LEN_TRIM(namr) == 0) THEN
         actio = 'NEW   '
       ELSE
         actio = 'RESTAR'
         error(5) = .FALSE.
       END IF
       error(8) = .FALSE.
     END IF

     !--------------------error checking
     !                    no datafile
     IF(error(2))WRITE(*,"(' ERROR: no input file name has been given')")
     !                    for restart
     IF(TRIM(actio) == 'RESTAR' .OR. TRIM(actio) == 'NSTRA0') THEN
       !                    -no restart file
       IF(error(5))WRITE(*,"(' ERROR: no restart file name has been given')")
       !                    reset error flag IF no restart
     ELSE
       error(5)=.FALSE.
     END IF

     !                    any error
     IF(ANY(error))THEN
       WRITE(*,1000)
       WRITE(55,1000)
       STOP
     END IF

   END IF
   IF (INDEX(output,'.@')/=0 .OR. INDEX(output,'.@')/=0) THEN
     WRITE(*,1001)
     WRITE(55,1001)
     STOP
   END IF

   IF (TRIM(actio)=='NEW') CALL del_old_files(1)   !Borra los ficheros de ejecuciones anteriores

   CALL openfi(55)      ! open report file
   IF (TRIM(actio) == 'NEW') THEN    !no difference presently
     CALL screin('BEGINP',cpui)
   ELSE
     CALL screin('RESTAR',cpui)
   END IF

   data_file = TRIM(prognm)//'.dat'
   CALL openfi(5)      ! 5 is the flag to open .DAT & .RSN files

   IF (TRIM(actio) == 'NEW') THEN
     CALL wrtrpt('BEGINP',narg=narg,params=params)
   ELSE
     CALL wrtrpt('RESTAR',narg=narg,params=params)
   END IF

   IF (actchk) THEN
     WRITE(*,"('*** CHECKING DATA INPUT FILE FOR FATAL ERRORS. PLEASE WAIT ...',/)")
     WRITE(55,"('*** CHECKING DATA INPUT FILE FOR FATAL ERRORS. PLEASE WAIT ...',/)")
   ELSE
     WRITE(*,"('*** PLEASE WAIT FOR THE PRELIMINARY CALCULATIONS TO BE COMPLETED ...',/)")
     WRITE(55,"('*** PLEASE WAIT FOR THE PRELIMINARY CALCULATIONS TO BE COMPLETED ...',/)")
   END IF
    memo = MIN(memo,max_memo)
 RETURN

  1000 FORMAT(' '/                                                                   &
      &' The PROGRAM should be executed with the following statement:'/              &
      &' exe_name [param] [INP=]file1 [OUT=gname] [FRQ=num1] [RSF=file2]'            &
      &/'           [MEM=num2] [FRS=num3] CHK '/                                     &
      &' * param : Is either "NEW" or "RES" or "NST" or "LST". (default=NEW)'/       &
      &' * file1 : Is the DATA input filename (no default).'/                        &
      &'           It should be composed of <GenericName>.<ext>'/                    &
      &' * file2 : Is the restart filename (necessary only IF "RES" or "NST"'/       &
      &'           is used). It should be composed of <GenericName>.<ext>'/          &
      &'           If <ext> is unknown the extension "rsf" can be used.'/            &
      &' * gname : Is the filename for output. (default=generic_name)'/              &
      &' * num1  : Is an INTEGER number for restart dump frequency in minutes.'/     &
      &'           (default=30)(Negative value for number of time steps frequency)'/ &
      &' * num2  : Is an INTEGER number for MEMory option (default=4).'/             &
      &' * CHK   : If this word is present data file is checked.'/                   &
      &' Parameters may be in any order and'/                              &
      &' key words can be written either in lower or uppercase letters.'/  &
      &' The following shorthand version may also be used:'/               &
      &' exec_name  file1'/                                                &
      &'  where all the above defaults apply')
 1001 FORMAT(' '/ &
      &' OUTput base filename must not include substring  .@  please modify')
 END SUBROUTINE beginp
