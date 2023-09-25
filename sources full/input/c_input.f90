 MODULE c_input
 ! data base and routine to manage input
 USE param_db,ONLY: midn, miln, mnam, mstl, mchr, mich, mlen, mttl
 USE lispa0
 !*** variables for listening and extracting values
 IMPLICIT NONE


   LOGICAL :: backs = .FALSE.        !TRUE to keep WORDS & PARAMS from previous call to LISTEN

   ! standart variables for reading reals RDFRRE

   INTEGER, PARAMETER :: maxnr = 40  !maximum number of REALS than can be read in a line
   INTEGER (kind=4) :: nreal         !actual number of REALs read
   REAL    (kind=8) :: reard(maxnr)  !values read

   ! standart variables for reading integers  RDFRIN

   INTEGER, PARAMETER :: maxni = 40   !maximum number of INTEGERS than can be read in a line
   INTEGER (kind=4)   :: nintg,   &   !actual number of INTEGERS read
                         intrd(maxni) !values read

   LOGICAL :: echo = .FALSE. ! DEFAULT echo .FALSE. to change it use: echo .TRUE.
   INTEGER (kind=4)   ::  ierrd       !flag indicating erro

 CONTAINS

   SUBROUTINE upcase(string)
   !*************************************************************************
   !
   !****this routine modifies a string, converting  a..z ==> A..Z
   !
   !*************************************************************************
   IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*),INTENT(INOUT):: string

   !--- Local variables
   INTEGER(kind=4):: i, k, n

   n = LEN_TRIM(string)
   DO i=1,n
     k = ICHAR(string(i:i))
     IF(k >= 97 .AND. k <= 122)string(i:i)=char(k-32)
   END DO

   RETURN
   END SUBROUTINE upcase

   SUBROUTINE listen(subr)
   !***********************************************************************
   !
   !***reads a string and interprets it as words and parameters.
   !
   !   - maximum number of words and parameters = maxwp.
   !   - only the first six characters of each word are decoded.
   !   - the underline characters '_' are discarted.
   !   - lower case letters are converted to upper case.
   !   - each word or parameter must be separated by ' ', '=', ':' or ','
   !   - a "comment" begins with '$', '!', '/'  or '\' in any place.
   !   - a line that has a comment beginning with '/' or '\'
   !     continuates in the next line.
   !   - words between "[" and "]" are units and are not decoded.
   !   - a line beginning with title is not decoded. title directive.
   !   - a line beginning with echo turns 'on' or 'off' the echo.
   !
   !***********************************************************************
   IMPLICIT NONE

     !***  arguments
     CHARACTER(len=*),INTENT(IN)   :: subr  !calling subroutine
     !--- Local variables
     REAL    (kind=8) digit
     INTEGER (kind=4) first,firsp,last,lastp,ptrwo,npptr,nwptr,leng,flag,nname
     CHARACTER :: letra, letra1
     LOGICAL newln,fastr

     !***  begin
     IF( backs )THEN    !if data exists from the previous call to listen
       backs = .FALSE.  !modify flag
       RETURN           !end of the task
     END IF

     nnwor = 0         ! initialize number of WORDS
     nname = 0         ! initialize number of NAMES
     nwptr = 0         ! initialize pointer to WORDS
     npptr = 0         ! initialize pointer to PARAMS
     ierrd = 0         ! error indicator
     words = ''        ! initializes WORDS
     param = 0.0d0     ! initializes PARAMETERS
     newln = .FALSE.   ! initializes continuation line to FALSE
     firsp = 0         ! initializes to zero just to have some initial value

     ! ***** use fast read first ******
     fastr = .TRUE.      !initializes
     CALL rdfrre('LISTEN',param,nnpar,maxnr,fastr)
     IF( fastr )THEN  !check is read was successful
       nwopa = nnpar  !number of WORDS/PARAMS = number of PARAMS
       RETURN         !end of the task
     END IF
     ! ***** WORDS and PARAMS are mixed
     nnpar = 0         ! initialize number of PARAMETERS
     param = 0.0d0     ! initializes PARAMETERS

     DO
       IF(((nnwor /= 0).OR.(nnpar /= 0))    &  ! Don't RETURN without answer
                      .AND. .NOT.newln) EXIT   ! CONTINUE reading IF / or \

       last  = 0                       ! initializes last processed character
       lastp = 0

       READ(ludat,"(A)",end=2,err=1) card     ! READ a card
       newln = .FALSE.                        ! initialize
       leng = LEN_TRIM(card)                  ! compute the length.

       outer : DO                             ! decode all the card.
         IF(last >= leng)EXIT                 ! all characters processed
         !  block to find the FIRST character of the Key-word
         first = last+1                       ! first non-processed character
         DO                          ! Skip over auxiliar characters
           IF(first > leng) EXIT outer        ! end of card reached EXIT loop
           letra = card(first:first)          ! present character
           !       skip auxiliar characters
           IF(letra /= '_' .AND.  &           ! Jump null CHARACTER (_)
              letra /= CHAR(9).AND.  &        ! jump tabulator
              letra /= ' ' .AND.  &           ! Jump separators ( =:,)
              letra /= '=' .AND.  &
              letra /= ':' .AND.  &
              letra /= ',' .AND.  &
              letra /= '['    ) EXIT
           !       skip over Units
           IF(letra  == '[') THEN  ! Jump Units (Ejem: [KN/m^2])
             DO
               first = first+1                ! skip previous character
               IF(first > leng) CALL runend & ! non-paired [], ERROR
                                ('LISTEN: UNITS MUST END WITH "]"    ')
               letra = card(first:first)      ! present character
               IF( letra == ']')EXIT          ! closing bracket found, EXIT loop
             END DO
           END IF
           first = first+1                    ! new first non-processed char
         END DO
         !         read key word
         IF(last == 0) firsp = first          ! SAVE first to PRINT card
         !  loop to find the LAST character of the Key-word
         last = first                         !initializes
         ! check if a label is input
         IF( letra == "'" .OR. letra == '"' )THEN
           IF( nnwor < 1)THEN  !if it begins with a label
             nnwor = 1
             nnpar = 1
           END IF
           letra1 = letra                       ! keep opening char
           DO                                   ! search for closing char
             last = last+1                      ! point to next character
             IF( last > leng )CALL runend('LISTEN: non paired label chars     ')
             letra = card(last:last)            ! present character
             IF( letra == letra1 )EXIT          ! closing char found
           END DO
           IF( last - first < 3)CALL runend('LISTEN: a null label is not valid  ')
           nname = nname + 1                    ! increase number of names
           names(nname) = ''                    ! initializes name blank
           names(nname) = card(first+1:MIN(last-1,mnam+first+1))
           param(nnwor) = - nname
         ELSE
           DO  ! search for the last character, and then exit
             letra = card(last:last)            ! present character
             IF( letra == ' ' .OR. &            ! Look for separator ( =:,).
                 letra == CHAR(9).OR.  &        ! jump tabulator
                 letra == '=' .OR. &             !standard separators
                 letra == ':' .OR. &             !
                 letra == ',' .OR. &             !
                 letra == '[' .OR. &             ! Look for units
                 letra == '$' .OR. &             ! Look for coment ($!).
                 letra == '!' .OR. &
                 letra == '/' .OR. &             ! Look for continuation (/\).
                 letra == '\' ) EXIT
             last = last+1                      ! point to next character
             IF( last > leng )EXIT              ! look for END of line
           END DO

           last = last - 1                      ! correct last position
           IF( letra == '$' .OR. letra == '!')  leng = last  !comments
           IF(letra == '/' .OR. letra == '\') THEN     ! continuation line
             leng  = last                        ! deal with continuation
             newln = .TRUE.                      ! set new line flag.
           END IF
           IF(last < first) EXIT               ! nothing read, EXIT loop
           lastp = last                        ! save last to print card
           READ(card(first:last),'(G20.0)',IOSTAT = flag) digit  !read a number
           IF (card(first:first)=='d' .OR.  &     !there is compiler error
               card(first:first)=='D' .OR.  &     !if the character strings
               card(first:first)=='e' .OR.  &     !begins with 'd' or 'e'
               card(first:first)=='E' ) flag=64   !and the rest is a number
           IF(flag == 0) THEN                  ! O.K. it is a parameter.
             nnpar = nnpar+1                   ! increase # of parameters
             npptr = npptr+1                   ! pointer to next parameter
             IF(npptr > maxwp) CALL runend  &  ! check Maximum allowed
                               ('LISTEN: TOO MANY PARAM IN COMMAND  ')
             IF(nwptr > npptr) npptr = nwptr   ! syncronizes words & params
             param(npptr) = digit              ! store param
           ELSE                                ! it is a word.
             nnwor = nnwor+1                   ! # of words
             nwptr = nwptr+1                   ! pointer to next word
             IF(nwptr > maxwp) CALL runend &   ! check Maximum allowed
                               ('LISTEN: TOO MANY WORDS IN COMMAND  ')
             IF(npptr >=  nwptr) nwptr=npptr+1 ! syncronizes words & params
             ptrwo = 0                         ! initializes word length
             DO                                ! loop to extract the WORD
               ptrwo = ptrwo+1                 ! updates word length
               ! if word processed or word length greater than 6, Exit loop
 !              IF (first > last .OR. ptrwo == midn+1) EXIT
               IF (first > last .OR. ptrwo == mstl+1) EXIT
               words(nwptr)(ptrwo:ptrwo) = card(first:first) ! store character
               first = first+1                 ! process next character
               DO   ! loop to Jump over Underline CHARACTERs
                 IF(card(first:first) /= '_') EXIT
                 first = first+1
               END DO
             END DO
             CALL upcase(words(nwptr))       ! convert to upper case.
           END IF

           IF(TRIM(words(1)) == 'TITLE' ) THEN        ! deal with title
             IF(echo .AND. (TRIM(subr) /= 'NOECHO')) &
               WRITE(lures,"(1X,A6,' <-- ',A)",ERR=9999) 'LISTEN',card(firsp:leng)
             first = last+2
             DO
               IF(card(first:first) /= ' ')EXIT  ! Remove trailing blank spaces
               first = first + 1
             END DO
             IF(leng < first) CALL runend ('LISTEN: BLANK IS ILEGAL HEADER     ')
             card = card(first:leng)              ! remove words(1) from card
             leng = leng-first+1
             firsp = 1
             lastp = leng
             param(1) = REAL(leng,8)              ! store length to print title
             EXIT outer                           ! break the DO
           END IF
         END IF
       END DO outer

       IF((TRIM(words(1)) == 'ECHO') .AND. (nnpar == 0) .AND. (nnwor == 2)) THEN  ! deal with echo
         IF(TRIM(words(2)) == 'OFF') THEN
           echo = .FALSE.
           IF(TRIM(subr) /= 'NOECHO')WRITE(lures,"(1X,A6,' <-- ',A)",ERR=9999)'LISTEN','ECHO OFF'
         ELSE
           echo = .TRUE.
           IF(TRIM(subr) /= 'NOECHO')WRITE(lures,"(1X,A6,' <-- ',A)",ERR=9999)'LISTEN','ECHO ON'
         END IF
         nnwor = 0                               ! forget all
         nwptr = 0
         words = ''
       ELSE                                      ! PRINT card
         IF((echo).AND.(firsp <= lastp) .AND. (TRIM(subr) /= 'NOECHO')) THEN
           IF(newln) THEN
             lastp = lastp+2
             card(lastp-1:lastp) = ' \\'
           END IF
           WRITE(lures,"(1X,A,' <-- ',A)",ERR=9999) TRIM(subr),card(firsp:lastp)
         END IF
       END IF

     END DO ! (nnwor /= 0) .OR. (nnpar /= 0) .AND. newln)

     nwopa  =  MAX(npptr,nwptr)

   RETURN
   !***  error:
  1 CALL runend('LISTEN: ERROR DETECTED WHEN READING')
  2 CALL runend('LISTEN: END OF FILE DETECTED       ')
  9999 CALL runen2('')
   END SUBROUTINE listen

   FUNCTION exists(fword,posit,value)
   !**********************************************************************
   !
   !****this routine verifys IF FWORD exists in WORDS
   !
   !**********************************************************************
   IMPLICIT NONE

     !***  arguments
     LOGICAL   exists
     CHARACTER(len=*),INTENT(IN) :: fword              !maximum length is 'midn'
     ! both optional arguments are unmodified if FWORD not found
     INTEGER (KIND=4), OPTIONAL, INTENT(OUT) :: posit  !possition in array WORDS
     REAL (KIND=8), OPTIONAL, INTENT(OUT) :: value     !associated value (PARAM)

     !     local
     INTEGER (kind=4)  iword,l

     !     begin

     exists = .FALSE.  ! initialize
     iword  = 0                                        !initializes pointer
     l = MIN(LEN_TRIM(fword),midn)                     !length of string to search

     DO
       IF( iword >= nwopa )EXIT
       iword  = iword + 1
       exists = (words(iword)(1:l) == fword(1:l))
       IF( exists ) THEN
         IF( PRESENT(posit)) posit = iword           !position in array
         IF( PRESENT(value)) value = param(iword)    !associated parameter
         EXIT
       END IF
     END DO

   RETURN
   END FUNCTION exists

   FUNCTION getint(fword,defau,texts)
   !**********************************************************************
   !
   !****this routine gets the INTEGER value associated with fword
   !
   !      - TEXTS(1:1) == '!' => COMPULSORY PARAMETER
   !
   !**********************************************************************
   IMPLICIT NONE

     INTEGER (kind=4)              :: getint
     !
     !***  arguments
     !
     CHARACTER(len=*), INTENT(IN)  :: fword
     INTEGER (kind=4), INTENT(IN)  :: defau
     CHARACTER(len=*), INTENT(IN) :: texts

     !     local
     INTEGER (kind=4)  iword, & !counter
                       i,j,l    !auxiliar
     LOGICAL    found
     CHARACTER(len=40) :: points='........................................'
     CHARACTER(len=7) :: spaces='       '

     !***  begin

     l = LEN_TRIM(texts)
     found = exists(fword,iword)
     IF(found) THEN
       getint = INT(param(iword))

     ELSE IF( texts(1:1) == '!')  THEN
       WRITE(lures,1,ERR=9999) TRIM(fword),texts(2:l)
       WRITE(*    ,1,ERR=9999) TRIM(fword),texts(2:l)
       CALL runend('GETINT: COMPULSORY PARAM. NOT FOUND')

     ELSE
       getint = defau

     END IF

     IF(texts(2:2) /= '-')THEN
       i = 40-l
       j = 7-LEN_TRIM(fword)
       WRITE (lures,"(9X,A,A,' ',A,A,' = ',i12)",ERR=9999) texts(2:l),points(1:i),TRIM(fword),spaces(1:j),getint
     END IF

   RETURN
   !***  format for error messagge
  1   FORMAT(/,4X,'*** ERROR: ',A,' = ',A,/,15X,'IS A COMPULSORY PARAMETER. SPECIFY IT !')
  9999 CALL runen2('')
   END FUNCTION getint

   FUNCTION getrea(fword,defau,texts)
   !**********************************************************************
   !
   !****this routine gets the REAL value associated with fword
   !
   !      - TEXTS(1:1) == '!' => COMPULSORY PARAMETER
   !
   !**********************************************************************
   IMPLICIT NONE

     REAL (kind=8) getrea

     !***  argument
     CHARACTER(len=*), INTENT(IN) :: fword
     REAL (kind=8),    INTENT(IN) :: defau
     CHARACTER(len=*), INTENT(IN) :: texts

     !     local
     INTEGER (kind=4)  iword, & !counter
                       i,j,l    !auxiliar
     LOGICAL    found
     CHARACTER(len=40) :: points='........................................'
     CHARACTER(len=7) :: spaces='       '

     !***  begin

     l = LEN_TRIM(texts)
     found = exists(fword,iword)
     IF(found) THEN
       getrea = param(iword)

     ELSE IF( texts(1:1) == '!')  THEN
       WRITE(lures,1,ERR=9999) TRIM(fword),texts(2:l)
       WRITE(*,    1,ERR=9999) TRIM(fword),texts(2:l)
       CALL runend('GETINT: COMPULSORY PARAM. NOT FOUND')

     ELSE

       IF( ABS(defau) > 1D-19 )THEN
         getrea = defau
       ELSE
         getrea = 0d0
       END IF

     END IF

     IF(texts(2:2) /= '-')THEN
       i = 40-l
       j = 7-LEN_TRIM(fword)
       WRITE (lures,"(9X,A,A,' ',A,A,' = ',G14.5)",ERR=9999) texts(2:l),points(1:i),TRIM(fword),spaces(1:j),getrea
     END IF

   RETURN
   !***  format for error messagge
  1   FORMAT(/,4X,'*** ERROR: ',A,' = ',A,/,15X,'IS A COMPULSORY PARAMETER. SPECIFY IT !')
  9999 CALL runen2('')
   END FUNCTION getrea

   SUBROUTINE rdfrin(subna,intrd,nintg,maxni)
   !***********************************************************************
   !
   !***reads a string and interprets it as integers (up to MAXNI per line)
   !
   !     It is assummed that only integers are read on the line !!!
   !***********************************************************************
   IMPLICIT NONE

     !--- Dummy variables
     INTEGER (kind=4),INTENT(IN)  :: maxni
     INTEGER (kind=4),INTENT(OUT) :: nintg, intrd(maxni)
     CHARACTER(len=*),INTENT(IN)  :: subna
     !--- Local varibles
     INTEGER (kind=4) i,last,leng,first,flag
     LOGICAL  newln

     !***  begin
     IF( backs )THEN  !is this possible? I mean a call to rdfrin after a BackSpace (FF)
       nintg = nnpar
       intrd(1:nintg) = INT(param(1:nnpar))
       backs = .FALSE.
       RETURN
     END IF
     intrd = 0
     nintg = 0

     DO
       READ(ludat,"(A)",end=2,err=1) card            ! READ a card
       leng = LEN_TRIM(card)                         ! calculate the length.
       IF(leng == 0 )CYCLE
       ! first search for comment sign $ or ! and eliminate remaining characters
       i = 1
       remark : DO
         IF( card(i:i) == '$' .OR. card(i:i) == '!' ) THEN
           leng = i
           DO
             leng = leng - 1
             IF(leng == 0) EXIT remark
             IF( card(leng:leng) /= ' ' ) EXIT remark
           END DO
         END IF
         IF( i == leng)EXIT
         i = i+1
       END DO remark

       IF( leng == 0 )CYCLE

       ! determines if it is a continuation line
       IF( card(leng:leng) == '/' .OR. card(leng:leng) == '\' ) THEN
         newln=.TRUE.
         DO
           leng = leng - 1
           IF(leng == 0) EXIT
           IF( card(leng:leng) /= ' ' ) EXIT
         END DO
       ELSE
         newln = .FALSE.
       END IF

       i = 0
       DO
         i = i + 1
         IF (i > leng) EXIT
         IF (card(i:i) /= ' ') THEN
           nintg = nintg+1
           first = i
           last  = leng
           DO
             i = i+1
             IF( i > leng) EXIT
             IF (card(i:i) == ' ') THEN
               i = i-1
               last = i
               EXIT
             END IF
           END DO

           IF(nintg <= maxni)THEN
             READ(card(first:last),'(I20)',IOSTAT=flag) intrd(nintg)
             IF (flag /= 0) THEN
               WRITE (lures,"(A)",ERR=9999) TRIM(subna)
               WRITE (lures,"(A)",ERR=9999) TRIM(card)
               CALL runend('RDFRIN: AN INTEGER WAS EXPECTED   ')
             END IF
           END IF
         END IF
       END DO
       IF (nintg > maxni )THEN
         WRITE(lures,*,ERR=9999)'  Warning !!! More INTEGERS than Expected.'
         WRITE(lures,*,ERR=9999) TRIM(card)
       END IF
       IF ( newln ) CYCLE
       EXIT
     END DO

   RETURN
   !     ***errors:
  1 CALL runend('RDFRIN: ERROR DETECTED WHEN READING')
  2 CALL runend('RDFRIN: END OF FILE DETECTED       ')
  9999 CALL runen2('')
   END SUBROUTINE rdfrin

   SUBROUTINE rdfrre(subna,reard,nreal,maxnr,fastr)
   !***********************************************************************
   !
   !***READS A STRING AND INTERPRETS IT AS REALS (UP TO MAXNR PER LINE)
   !
   !     It is assummed that only reals are read on the line !!!
   !***********************************************************************
   IMPLICIT NONE

     !--- Dummy variables
     INTEGER (kind=4),INTENT(IN)  :: maxnr
     INTEGER (kind=4),INTENT(OUT) :: nreal
     REAL    (kind=8),INTENT(OUT) :: reard(maxnr)
     CHARACTER(len=*),INTENT(IN)  :: subna
     LOGICAL, INTENT(IN OUT) :: fastr    !fast read
     !--- Local variables
     INTEGER (kind=4) i,last,leng,first,flag,nrdln
     LOGICAL  newln

     !***  begin

     reard = 0  !initializes real parameters read
     nreal = 0  !initializes number of parameters read
     nrdln = 0  !initializes number of phisycal lines read

     DO
       READ(ludat,"(A)",end=2,err=1) card            ! READ a card
       leng = LEN_TRIM(card)                         ! calculate the length.
       nrdln = nrdln + 1
       IF(leng == 0 )CYCLE
       ! first search for comment sign $ and eliminate remaining characters
       i = 1
       remark : DO
         IF( card(i:i) == '$' ) THEN
           leng = i
           DO
             leng = leng - 1
             IF(leng == 0) EXIT remark
             IF( card(leng:leng) /= ' ' ) EXIT remark
           END DO
         END IF
         IF( i == leng)EXIT
         i = i+1
       END DO remark

       IF( leng == 0 )CYCLE

       ! determines if it is a continuation line
       IF( card(leng:leng) == '/' .OR. card(leng:leng) == '\' ) THEN
         newln=.TRUE.
         DO
           leng = leng - 1
           IF(leng == 0) EXIT
           IF( card(leng:leng) /= ' ' ) EXIT
         END DO
       ELSE
         newln = .FALSE.
       END IF

       i = 0
       DO
         i = i + 1
         IF (i > leng) EXIT
         IF (card(i:i) /= ' ') THEN
           nreal = nreal+1
           first = i
           last  = leng
           DO
             i = i+1
             IF( i > leng) EXIT
             IF (card(i:i) == ' ') THEN
               i = i-1
               last = i
               EXIT
             END IF
           END DO

           IF(nreal <= maxnr)THEN
             READ(card(first:last),'(f20.0)',IOSTAT=flag) reard(nreal)
             IF (flag /= 0) THEN
               IF( fastr ) THEN
                 fastr = .FALSE.
                 DO i=1,nrdln
                   BACKSPACE (ludat)
                 END DO
                 nwopa = nnpar
                 RETURN
               ELSE
                 WRITE (lures,"(A)",ERR=9999) TRIM(subna)
                 WRITE (lures,"(A)",ERR=9999) TRIM(card)
                 CALL runend('RDFRRE: A REAL NUMBER WAS EXPECTED')
               END IF
             END IF
           END IF
         END IF
       END DO
       IF (nreal > maxnr ) THEN
         WRITE(lures,*,ERR=9999)'  Warning !!! More REALS than Expected.'
         WRITE(lures,*,ERR=9999) TRIM(card)
       END IF
       IF ( newln ) CYCLE
       EXIT
     END DO
   RETURN
   !     ***errors:
  1 CALL runend('RDFRRE: ERROR DETECTED WHEN READING')
  2 CALL runend('RDFRRE: END OF FILE DETECTED       ')
  9999 CALL runen2('')
   END SUBROUTINE rdfrre

   SUBROUTINE openfi(nfile,rest,narch,flag,oflag)
   !*********************************************************************
   !
   !     OPEN files to be used
   !
   !*********************************************************************

   USE ifport  !Intel Compiler
   USE name_db, ONLY: input, output, prognm, rsfprj, rsfstg, rsfnam, data_file,gen_data_file
   USE gvar_db, ONLY: actchk
   IMPLICIT NONE

     ! Dummy arguments
     INTEGER(kind=4),INTENT(IN):: nfile
     CHARACTER(len=*),OPTIONAL:: rest
     INTEGER(kind=4),INTENT(IN),OPTIONAL:: flag
     INTEGER(kind=4),INTENT(OUT),OPTIONAL:: narch, oflag
     ! Functions
     CHARACTER(len=mich):: inttoch
     ! Local variables
     INTEGER(kind=4),SAVE :: sfile  !21 -- 40
     CHARACTER(len=mstl):: auxi(13)
     CHARACTER(len=7):: st_key  ! 'UNKNOWN'  or 'SCRATCH'
     CHARACTER(len=6):: po_key  ! 'ASIS  '  or 'APPEND'

     INTEGER(kind=4):: i, j, fpatn, ist, dummy, lo, lr
     CHARACTER(len=mlen):: stval,fpath
     CHARACTER(len=mlen):: string
     LOGICAL:: endst

     !!varname='STAMP_DTS'
     !!fpatn = GETENVQQ(varname(1:9), STVAL)
     stval = '.'
     fpatn = 1
     fpath = stval(1:fpatn)//'\'
     fpatn = fpatn + 1

  ! Status key
  st_key = 'UNKNOWN'    !default
  IF( actchk ) st_key = 'SCRATCH' !for active checks
  ! Position key
  po_key = 'ASIS  '     !default
  IF(PRESENT(flag))THEN
    IF(flag == 1) po_key = 'APPEND'
  END IF

  SELECT CASE (nfile)
     CASE (1)    !       Data file for implicit FE program: springback
       OPEN( 1,FILE=fpath(1:fpatn)//TRIM(output)//'.rel',FORM='formatted',STATUS=st_key)

     CASE (2)    !       time information file
       OPEN( 2,FILE=fpath(1:fpatn)//TRIM(output)//'.tim',FORM='formatted',STATUS=st_key)

     CASE (5)    !used as flag to open ASCII files
        !(3) (.rsn) : output may be full or not
        !(5) (.dat) : generated input data file
        !(58) (.deb) : debug file for developers

      IF( gen_data_file )THEN
        OPEN(7,FILE=fpath(1:fpatn)//TRIM(input),FORM='formatted',STATUS='unknown')   !data input file
        OPEN(3,FILE=fpath(1:fpatn)//TRIM(data_file),FORM='formatted',STATUS='old',  &
                                  IOSTAT=ist)                                        !check if PROGRAM.dat file exist
        IF (ist==0) CLOSE(3,STATUS='delete')                                         !and delete
        OPEN(3,FILE=fpath(1:fpatn)//TRIM(data_file),FORM='formatted',STATUS='new')   !generated file PROGRAM.dat
        CALL genfil (7,fpatn,fpath)
        CLOSE (3)
        CLOSE (7)
        OPEN(ludat,FILE=fpath(1:fpatn)//TRIM(data_file),FORM='formatted',STATUS='old')
      ELSE
        OPEN(ludat,FILE=fpath(1:fpatn)//TRIM(input),FORM='formatted',STATUS='old')
      END IF
      OPEN(lures,FILE=fpath(1:fpatn)//TRIM(output)//'.rsn',FORM='formatted',STATUS='unknown')

      OPEN(58,FILE=fpath(1:fpatn)//TRIM(output)//'.deb',FORM='formatted',STATUS='unknown')

     CASE (6)   !to redirect screen messages to file
        OPEN( 6,FILE=fpath(1:fpatn)//TRIM(output)//'.scr',FORM='formatted',            &
                STATUS='unknown')

     !Acceleration files
     CASE (7)
        OPEN( 7,FILE=fpath(1:fpatn)//'accel.x',FORM='FORMATTED',STATUS='OLD')
     CASE (8)
        OPEN( 8,FILE=fpath(1:fpatn)//'accel.y',FORM='FORMATTED',STATUS='OLD')
     CASE (9)
        OPEN( 9,FILE=fpath(1:fpatn)//'accel.z',FORM='FORMATTED',STATUS='OLD')

     CASE (10)  !Auxiliar ASCII file
       OPEN(10,FILE=fpath(1:fpatn)//TRIM(rest),FORM='formatted',STATUS=st_key)

     CASE (11:15,19:20,70)   !           files for plots of displ of internal or contact forces
                             !           velocities and accelerations

      OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.p'//TRIM(inttoch(nfile,2)),  &
                 FORM='unformatted',STATUS=st_key,ACCESS='sequential',POSITION=po_key,SHARED)

     CASE (16:18)    !Files for GLOBAL OUTPUT
       sfile = 20 !necessary to initialize again for a new strategy
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.f'//TRIM(inttoch(nfile,2)), &
            FORM='unformatted',STATUS=st_key,POSITION=po_key,SHARED)

     CASE (21)  !                 files for plots of stresses
       sfile = sfile+1
       IF( sfile > 40 )THEN
         WRITE(lures,*,ERR=9999) 'TOO many sets with required stresses Max=20'
         CALL runend('TOO many sets with required stresses Max=20')
       END IF
       narch = sfile
       OPEN(narch,FILE=fpath(1:fpatn)//TRIM(output)//'.p'//TRIM(inttoch(nfile,2)),  & !  nfile->narch?
            FORM='unformatted',STATUS=st_key,POSITION=po_key,SHARED)

     CASE (22:40)  !                 files for plots of stresses
       WRITE(lures,*,ERR=9999) 'OPENFI: tries to open file reserved for stresses',nfile
       CALL runend('OPENFI: internal error with file   ')

     CASE (41)    !               files for plots of total forces over contact surfaces
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.ctc',FORM='unformatted',      &
             STATUS=st_key,POSITION=po_key,SHARED)


     CASE (42)  !Unit 42 - saving drawbead force history
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.p42',FORM='unformatted',      &
            STATUS=st_key,POSITION=po_key,SHARED)

     CASE (43)  !Unit 43 - output of drawbead lines
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.dbl', FORM='unformatted',     &
            STATUS=st_key,POSITION=po_key,SHARED)

     CASE (44)  !Nodal contact forces and distances on slave surface
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.fc4',FORM='unformatted',      &
            STATUS=st_key,POSITION=po_key,SHARED)

     CASE (45)  !Nodal friction work on master surface
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.fw4',FORM='unformatted',      &
            STATUS=st_key,POSITION=po_key,SHARED)

     CASE (46)  ! Binary export file
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(rest),FORM='unformatted',STATUS=st_key)

     CASE (47:49) !Any auxiliar file
       OPEN(nfile,FORM='unformatted',STATUS='scratch')

     CASE (50)  ! Dumping Files for restart
       rsfnam = TRIM(output)//TRIM(rest)
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(rsfnam),FORM='unformatted',STATUS=st_key)

       IF (flag == 0) THEN  !Restart file at the end of the strategy
         WRITE(50,ERR=9999) .TRUE.

       ELSE IF (flag == 1) THEN  !Restart file at the middle of the strategy
         WRITE(50,ERR=9999) .FALSE.

         !--- Copy output file for restart process
         oflag = 0
         !Copy outpuf file *.f16 from restart
         dummy = copy_file( fpath(1:fpatn)//TRIM(output)//'.f16',            &
                            fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.f16' )
         IF (dummy /= 0) oflag=-1   !Verify if command could not be done

         !Copy outpuf file *.f17 from restart
         dummy = copy_file( fpath(1:fpatn)//TRIM(output)//'.f17',            &
                            fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.f17' )
         IF (dummy /= 0) oflag=-1   !Verify if command could not be done

         ! !copy outpuf file *.f18 from restart
         ! dummy = if_exist_copy_file( fpath(1:fpatn)//TRIM(output)//'.f18 '     &
         !                              fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.f18' )
         ! IF (dummy /= 0) oflag=-1  !Verify if command could not be done

         !copy outpuf file *.ctc from restart
         dummy = if_exist_copy_file( &
                    fpath(1:fpatn)//TRIM(output)//'.ctc ',  &
                    fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.ctc'  )
         IF (dummy /= 0) oflag=-1  !Verify if command could not be done

         dummy = if_exist_copy_file( &
                    fpath(1:fpatn)//TRIM(output)//'.fc4 ',  &
                    fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.fc4 '  )
         IF (dummy /= 0) oflag=-1  !Verify if command could not be done

         dummy = if_exist_copy_file( &
                    fpath(1:fpatn)//TRIM(output)//'.fw4 ',  &
                    fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.fw4 '  )
         IF (dummy /= 0) oflag=-1  !Verify if command could not be done


         DO i=11,15
            !copy outpuf file *.p?? from restart
            dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(output)//'.p'//TRIM(inttoch(i,2)), &
                 fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.p'//TRIM(inttoch(i,2)) )
            IF (dummy /= 0) oflag=-1  !Verify if command could not be done
         END DO
         DO i=19,40
            !copy outpuf file *.p?? from restart
            dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(output)//'.p'//TRIM(inttoch(i,2)), &
                 fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.p'//TRIM(inttoch(i,2)) )
            IF (dummy /= 0) oflag=-1  !Verify if command could not be done
         END DO
         ! temperatures
         dummy = if_exist_copy_file( &
              fpath(1:fpatn)//TRIM(output)//'.p70',       &
              fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.p70' )
         IF (dummy /= 0) oflag=-1  !Verify if command could not be done
         !
         dummy = if_exist_copy_file( &
              fpath(1:fpatn)//TRIM(output)//'.p42',       &
              fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.p42' )
         IF (dummy /= 0) oflag=-1  !Verify if command could not be done

         dummy = if_exist_copy_file( &
              fpath(1:fpatn)//TRIM(output)//'.dbl',       &
              fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'.dbl' )
         IF (dummy /= 0) oflag=-1  !Verify if command could not be done


       END IF

     CASE (51)  !Files for restart (read)

       string = TRIM(rest)
       CALL upcase(string)
       IF (INDEX(string,'.RSF') > 0) THEN
         OPEN(9999,FILE=fpath(1:fpatn)//TRIM(rest),FORM='FORMATTED',STATUS='OLD',IOSTAT=ist)
         IF (ist /= 0) CALL runen3('Unknown '''//TRIM(rest)//''' file.')
         READ(9999,"(A)",IOSTAT=ist) rest
         IF (ist /= 0) CALL runen3('Bad name in '''//TRIM(rest)//''' file.')
         CLOSE(9999,STATUS='keep')
       END IF

       !Open restart file
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(rest),FORM='unformatted',STATUS='old',IOSTAT=ist)
       IF (ist /= 0) CALL runen3('Unknown '''//TRIM(rest)//''' file.')
       READ(51) endst

       !--- Copy output file for restart process
       lo = LEN_TRIM(output)
       lr = lfile(rest)
       IF (endst) THEN
         IF (output(1:lo) /= rest(1:lr)) THEN
           lo = LEN_TRIM(rest) - 4   !Elimina la extension
           lr = lr + 1
           string = TRIM(output)//rest(lr:lo)
           dummy = copy_file( fpath(1:fpatn)//rest(1:lo)//'.f16',    &
                              fpath(1:fpatn)//TRIM(string)//'.f16'    )

           dummy = copy_file( fpath(1:fpatn)//rest(1:lo)//'.f17',    &
                              fpath(1:fpatn)//TRIM(string)//'.f17'    )

          !dummy = copy_file( fpath(1:fpatn)//rest(1:lo)//'.f18',    &
          !                 fpath(1:fpatn)//TRIM(string)//'.f18'    )

           dummy = if_exist_copy_file( &
                  fpath(1:fpatn)//rest(1:lo)//'.ctc ',  &
                  fpath(1:fpatn)//TRIM(string)//'.ctc'  )

           dummy = if_exist_copy_file( &
                  fpath(1:fpatn)//rest(1:lo)//'.fc4 ',   &
                  fpath(1:fpatn)//TRIM(string)//'.fc4'   )

           dummy = if_exist_copy_file( &
                  fpath(1:fpatn)//rest(1:lo)//'.fw4 ',   &
                  fpath(1:fpatn)//TRIM(string)//'.fw4'   )

           DO i=11,15
              dummy = if_exist_copy_file( &
                   fpath(1:fpatn)//rest(1:lo)//'.p'//TRIM(inttoch(i,2)),  &
                   fpath(1:fpatn)//TRIM(string)//'.p'//TRIM(inttoch(i,2)))
           END DO
           DO i=19,40
              dummy = if_exist_copy_file( &
                   fpath(1:fpatn)//rest(1:lo)//'.p'//TRIM(inttoch(i,2)),  &
                   fpath(1:fpatn)//TRIM(string)//'.p'//TRIM(inttoch(i,2)))
           END DO
           dummy = if_exist_copy_file( &
                fpath(1:fpatn)//rest(1:lo)//'.p70',   &
                fpath(1:fpatn)//TRIM(string)//'.p70')

           dummy = if_exist_copy_file( &
                fpath(1:fpatn)//rest(1:lo)//'.p42',   &
                fpath(1:fpatn)//TRIM(string)//'.p42')

           dummy = if_exist_copy_file( &
                fpath(1:fpatn)//rest(1:lo)//'.dbl',   &
                fpath(1:fpatn)//TRIM(string)//'.dbl')

         END IF
       ELSE
         IF (output(1:lo) /= rest(1:lr)) THEN
           lo = LEN_TRIM(rest) - 4   !Elimina la extension
           lr = lr + 1
           string = TRIM(output)//rest(lr:lo)
         ELSE
           lo = LEN_TRIM(rest) - 4   !Elimina la extension
           string = rest(1:lo)
         END IF

         !Copy outpuf file *.f16 from restart
         dummy = copy_file( fpath(1:fpatn)//TRIM(rest)//'.f16 ',         &
                            fpath(1:fpatn)//TRIM(string)//'.f16'         )

         !Copy outpuf file *.f17 from restart
         dummy = copy_file( fpath(1:fpatn)//TRIM(rest)//'.f17 ',         &
                            fpath(1:fpatn)//TRIM(string)//'.f17'         )

         !Copy outpuf file *.f18 from restart
         dummy = if_exist_copy_file( fpath(1:fpatn)//TRIM(rest)//'.f18 ',         &
                            fpath(1:fpatn)//TRIM(string)//'.f18'         )

         !Copy outpuf file *.ctc from restart
         dummy = if_exist_copy_file( fpath(1:fpatn)//TRIM(rest)//'.ctc ',         &
                            fpath(1:fpatn)//TRIM(string)//'.ctc'         )

         !Copy outpuf file *.fc4 from restart
         dummy = if_exist_copy_file( fpath(1:fpatn)//TRIM(rest)//'.fc4 ',         &
                            fpath(1:fpatn)//TRIM(string)//'.fc4'         )

         !Copy outpuf file *.fw4 from restart
         dummy = if_exist_copy_file( fpath(1:fpatn)//TRIM(rest)//'.fw4 ',         &
                            fpath(1:fpatn)//TRIM(string)//'.fw4'         )

         DO i=11,15
            dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(rest)//'.p'//TRIM(inttoch(i,2)), &
                 fpath(1:fpatn)//TRIM(string)//'.p'//TRIM(inttoch(i,2)))
         END DO
         DO i=19,40
            dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(rest)//'.p'//TRIM(inttoch(i,2)), &
                 fpath(1:fpatn)//TRIM(string)//'.p'//TRIM(inttoch(i,2)))
         END DO
         dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(rest)//'.p70',     &
                 fpath(1:fpatn)//TRIM(string)//'.p70')
         dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(rest)//'.p42',     &
                 fpath(1:fpatn)//TRIM(string)//'.p42')
         dummy = if_exist_copy_file( &
                 fpath(1:fpatn)//TRIM(rest)//'.dbl',     &
                 fpath(1:fpatn)//TRIM(string)//'.dbl')

       END IF

     CASE (52)  !Files for plots of volume and pressure of follower load
      OPEN(52,FILE=fpath(1:fpatn)//TRIM(output)//'.vol',FORM='unformatted',     &
            STATUS=st_key,POSITION=po_key,SHARED)

     !Batch file for GiD
     CASE (53)
       OPEN(53,FILE=fpath(1:fpatn)//TRIM(output)//'.bgi',FORM='formatted',       &
            STATUS=st_key,SHARED)

     CASE (54)   !New mesh data (after remeshing)
       OPEN(54,FILE=fpath(1:fpatn)//TRIM(output)//'.nms',FORM='formatted',STATUS='OLD')

     CASE (55)  !report file
       OPEN(55,FILE=fpath(1:fpatn)//TRIM(output)//'.rep',FORM='formatted',              &
            ACCESS='SEQUENTIAL',STATUS='unknown',POSITION='APPEND',SHARED)

     CASE (56)   !Background mesh for remeshing
       OPEN(56,file=fpath(1:fpatn)//'mesh.dat',FORM='formatted',STATUS='OLD')
       DO j=1,13
         READ(56,"(A)") auxi(j)
       END DO
       CLOSE(56)
       OPEN(56,file=fpath(1:fpatn)//'mesh.dat',FORM='formatted',STATUS='OLD')
       WRITE(56,"(A)",ERR=9999) TRIM(auxi(1))
       WRITE(56,*,ERR=9999) fpath(1:fpatn)//TRIM(output)//'.bgm'
       WRITE(56,*,ERR=9999) fpath(1:fpatn)//TRIM(output)//'.nms'
       DO j=4,13
         WRITE(56,"(A)",ERR=9999) TRIM(auxi(j))
       END DO
       CLOSE(56)
       OPEN(56,FILE=fpath(1:fpatn)//TRIM(output)//'.bgm',FORM='formatted',STATUS=st_key)

     CASE(57)   ! List of node penetrations
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(rest)//'.pai',FORM='formatted',               &
            ACCESS='SEQUENTIAL',POSITION='APPEND',STATUS='unknown',SHARED)

     CASE (62)  !--- Report file
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(output)//'.rpt',FORM='formatted',             &
            ACCESS='SEQUENTIAL',POSITION='APPEND',STATUS=st_key,SHARED)

     CASE (63)  !--- Last restart file
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(rest)//'.rsf',FORM='formatted',          &
            ACCESS='SEQUENTIAL',STATUS='OLD',SHARED,IOSTAT=ist)
       rest = ''
       IF (ist == 0) THEN
         READ(nfile,"(A)") rest
         CLOSE(nfile)
       END IF

     !--- import data file
     CASE (71:90)
       OPEN(nfile,FILE=fpath(1:fpatn)//TRIM(rest),FORM='unformatted',STATUS='old', &
            IOSTAT=ist)
       IF( PRESENT(oflag) )oflag = ist

     CASE DEFAULT
       OPEN(nfile,FILE=TRIM(output)//TRIM(rest))
       WRITE(lures,"(' unknown file ',i3,' OPEN ')",ERR=9999) nfile
     END SELECT

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE openfi


   SUBROUTINE genfil (nr,fpatn,fpath)
   !     join DATA files into 'PROGRAM'.dat
   IMPLICIT none

     !        routine argument
     INTEGER(kind = 4),INTENT(IN):: nr,fpatn
     CHARACTER(len=mlen),INTENT(IN):: fpath
     !        local variables
     CHARACTER (len= 1) :: c
     CHARACTER (len= 8) :: form
     CHARACTER (len=miln) ::lin
     INTEGER (kind=4) :: l

     CALL randwr(c,lin,form,l,nr,fpatn,fpath)

   RETURN
   END SUBROUTINE genfil

   RECURSIVE SUBROUTINE randwr(c,lin,form,l,nr,fpatn,fpath)
   !
   IMPLICIT none

     INTEGER (kind=4),INTENT(IN) ::  nr,fpatn
     CHARACTER(len=mlen),INTENT(IN) :: fpath
     CHARACTER(len= 1),INTENT(IN OUT) :: c
     CHARACTER(len=*),INTENT(IN OUT) :: lin
     CHARACTER(len=*),INTENT(IN OUT) :: form
     INTEGER (kind=4),INTENT(IN OUT) :: l
     !     new file
     INTEGER (kind=4) :: nr1

     DO
       READ(nr,'(a1,A)',end = 1)c,lin
       SELECT CASE (c)
       CASE ('@')
         l = SCAN(lin,' ') - 1
         nr1 = nr+1
         OPEN (nr1,FILE=fpath(1:fpatn)//lin(1:l),FORM='formatted',STATUS='old')
         CALL randwr(c,lin,form,l,nr1,fpatn,fpath)
         CLOSE (nr1)
       CASE DEFAULT
         l = MAX(LEN_TRIM(lin),1)
 !        form = '(a1,a'//TRIM(inttoch(l,2))//')'
 !        WRITE(3,form,ERR=9999) c,lin(1:l)
         WRITE(3,"(A1,A)",ERR=9999) c,lin(1:l)
       END SELECT
     END DO

   1 RETURN
   9999 CALL runen2('')
   END SUBROUTINE randwr


   SUBROUTINE del_old_files(flag,rest)
   !*********************************************************************
   !
   !     DELETE old files to be used
   !
   !*********************************************************************
   USE ifport
   USE name_db, ONLY: output
   IMPLICIT NONE

     ! Dummy arguments
     INTEGER(kind=4),INTENT(IN) :: flag
     CHARACTER(len=*),OPTIONAL:: rest
     ! Local variables

     INTEGER(kind=4)::  dummy,fpatn
     CHARACTER(len=mlen):: varname
     CHARACTER(len=mlen):: stval,fpath

!#if UNIX
!     fpath = './'
!     fpatn = 2
!#else
     varname='STAMP_DTS'
     fpatn = GETENVQQ(varname(1:9), STVAL)
     stval = '.'
     fpatn = 1
     fpath = stval(1:fpatn)//'\'
     fpatn = fpatn + 1
!#endif


       SELECT CASE (flag)
       CASE (0)   !delete post-process files (GiD) to be newly generated
         ! CALLED from OUTDYN
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'*.post.*', &
                                      fpath(1:fpatn)//TRIM(output)//'*.post.*'  )
       CASE (1)   !delete all old files for a new problem
         ! CALLED from BEGINP
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.rep', &
                                      fpath(1:fpatn)//TRIM(output)//'.rep'   )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.rsn', &
                                      fpath(1:fpatn)//TRIM(output)//'.rsn'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.deb', &
                                      fpath(1:fpatn)//TRIM(output)//'.deb'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.rsf', &
                                      fpath(1:fpatn)//TRIM(output)//'.rsf'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//'FLC_*', &
                                      fpath(1:fpatn)//'FLC_*' )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.@1.???', &
                                      fpath(1:fpatn)//TRIM(output)//'.@1.???'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.@1.???.*', &
                                      fpath(1:fpatn)//TRIM(output)//'.@1.???.*'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.@1_f.???.*', &
                                      fpath(1:fpatn)//TRIM(output)//'.@1_f.???.*'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.@1-PAIR-*.???', &
                                      fpath(1:fpatn)//TRIM(output)//'.@1-PAIR-*.???'  )

       CASE (2)   !delete old files for a new strategy
         ! CALLED from NEWSTR
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.???', &
                                      fpath(1:fpatn)//TRIM(output)//'.???'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'.???.*',&
                                      fpath(1:fpatn)//TRIM(output)//'.???.*'  )
         dummy = if_exist_remove_file(fpath(1:fpatn)//TRIM(output)//'-PAIR-*.???',&
                                      fpath(1:fpatn)//TRIM(output)//'-PAIR-*.???'  )

       CASE (3) !Delete old restart files
            dummy = if_exist_remove_file( &
                 fpath(1:fpatn)//TRIM(output)//TRIM(rest),    &
                 fpath(1:fpatn)//TRIM(output)//TRIM(rest)//'*' )

       CASE DEFAULT
         CALL runend('del_old_files: non valid task, programming error')
       END SELECT

   END SUBROUTINE del_old_files

   SUBROUTINE manage_last_file_names (flag)
   USE ifport
   USE name_db, ONLY: output, rsfprj, rsfstg, rsfnam, namr
   IMPLICIT NONE

     ! dummy variables
     INTEGER(kind=4) :: flag

     ! Local variables

     INTEGER(kind=4)::  fpatn
     CHARACTER(len=mlen):: stval,fpath
     CHARACTER(len=mlen):: varname

!#if UNIX
!     fpath = './'
!     fpatn = 2
!#else
     varname='STAMP_DTS'
     fpatn = GETENVQQ(varname(1:9), STVAL)
     stval = '.'
     fpatn = 1
     fpath = stval(1:fpatn)//'\'
     fpatn = fpatn + 1
!#endif

     SELECT CASE (flag)
     CASE (1)
         !File with the name of the last strategy
         !Strategy level
         OPEN(9999,FILE=fpath(1:fpatn)//TRIM(rsfprj)//TRIM(rsfstg)//'.rsf',    &
              FORM='FORMATTED',STATUS='unknown')
         WRITE(9999,"(A)",ERR=9999) TRIM(rsfnam)
         CALL flushf(9999)
         CLOSE(9999,STATUS='keep')
         !Project level
         OPEN(9999,FILE=fpath(1:fpatn)//TRIM(rsfprj)//'.rsf',FORM='FORMATTED',STATUS='unknown')
         WRITE(9999,"(A)",ERR=9999) TRIM(rsfnam)
         CALL flushf(9999)
         CLOSE(9999,STATUS='keep')
     CASE (2)
         !File with the name of the last strategy
         OPEN(9999,FILE=fpath(1:fpatn)//TRIM(output)//'.rsf',FORM='FORMATTED',STATUS='unknown')
         WRITE(9999,"(A)",ERR=9999) TRIM(namr)
         CLOSE(9999,STATUS='keep')
         CALL flushf(9999)
     END SELECT
     RETURN
     9999 CALL runen2('')
   END SUBROUTINE manage_last_file_names

   SUBROUTINE rdreqs ( ngaus, nreqs, ngrqs, iwrit)
   !******************************************************************
   !
   !*** READ requested gauss points for output
   !
   !******************************************************************
   IMPLICIT NONE

     !--- Dummy variables
     INTEGER (kind=4), INTENT(IN) ::  ngaus, iwrit, nreqs
     INTEGER (kind=4), POINTER :: ngrqs(:)
     !--- Local variables
     INTEGER (kind=4) ::  i

     IF( nreqs <= 0 ) RETURN
     IF( ASSOCIATED(ngrqs) )DEALLOCATE(ngrqs)
     ALLOCATE(ngrqs(nreqs))
     CALL listen('RDREQS')

     IF (exists('GAUSSP')) THEN
       ! IF(nreqs == 0)CALL runend('RDREQS: NREQS EQUAL TO 0           ')
       IF(iwrit == 1) WRITE(lures,"(//'  REQUESTED GAUSS POINTS FOR OUTPUT', &
                                    /,'    ELEMENT      POINT '/)",ERR=9999)
       CALL rdfrin('RDREQS',ngrqs,nintg,nreqs)

       IF(iwrit == 1) WRITE(lures,"(2i10/)",ERR=9999) (((ngrqs(i)-1)/ngaus+1), &
                                  MOD(ngrqs(i)-1,ngaus)+1,i=1,nreqs)
     ELSE
       CALL runend('RDREQS: GAUSS_POINTS card expected ')
     END IF

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE rdreqs

   FUNCTION getchar(fword,defau,texts)
   ! this routine gets the CHARACTER (or INTEGER) value associated with fword
   !   - INTEGER is converted to CHARACTER
   !   - presently INTEGER < 100 (due to use of ch)
   !   - TEXTS(1:1) == '!' => COMPULSORY PARAMETER
   IMPLICIT NONE

     CHARACTER(len=mchr)           :: getchar

     ! arguments

     CHARACTER(len=*), INTENT(IN)  :: fword
     CHARACTER(len=*), INTENT(IN)  :: defau
     CHARACTER(len=*), INTENT(IN)  :: texts

     !--- Functions
     CHARACTER(len=mich):: inttoch

     ! local

     INTEGER (kind=4)   iword  !counter
     LOGICAL    found

     ! begin

     found = exists(fword,iword)
     IF (found) THEN
       IF (param(iword) == 0) THEN                         !if no associated number
         getchar = TRIM(words(iword+1)(1:MIN(mstl,mchr)))  !label is the next word
       ELSE                                                !label is a two digit number
         getchar = ''
         getchar = getchar(1:mchr-2)//TRIM(inttoch(INT(param(iword)),2)) !fill with blanks
       END IF

     ELSE IF( texts(1:1) == '!')  THEN
       WRITE(lures,1,ERR=9999) TRIM(fword),texts(2:LEN_TRIM(texts))
       WRITE(*    ,1,ERR=9999) TRIM(fword),texts(2:LEN_TRIM(texts))
       CALL runend('GETCHAR: COMPULSORY PARAM. NOT FOUND')

     ELSE
       getchar = defau

     END IF

     WRITE (lures,"(9X,A,'.... ',A,' = ',A)",ERR=9999) texts(2:LEN_TRIM(texts)),     &
            TRIM(fword), getchar

   RETURN
   ! format for error messagge
  1   FORMAT(/,4X,'*** ERROR: ',A,' = ',A,/,15X,'IS A COMPULSORY PARAMETER. SPECIFY IT !')
  9999 CALL runen2('')
   END FUNCTION getchar

   FUNCTION get_name(kword,founds,texts,defau,posin,posout,stype)
   ! this routine gets the NAME (label) value associated with key-word
   IMPLICIT NONE

     CHARACTER(len=mnam) :: get_name

     ! arguments
     CHARACTER(len=midn),INTENT(IN),OPTIONAL :: kword   !associated key-word
     LOGICAL, INTENT(OUT), OPTIONAL :: founds           !if key word is found
     CHARACTER(len=*), INTENT(IN), OPTIONAL  :: defau   !default value if key word not found
     CHARACTER(len=*), INTENT(IN), OPTIONAL  :: texts   !explanation text
     CHARACTER(len=*), INTENT(IN), OPTIONAL  :: stype   !set type to add to numbers
     INTEGER (kind=4), INTENT(IN), OPTIONAL  :: posin    !position in word list
     INTEGER (kind=4), INTENT(OUT), OPTIONAL :: posout  !position in word list
     !--- Functions
     CHARACTER(len=mich):: inttoch
     ! local
     INTEGER (kind=4)  iword, &  !counter
                       ival      !param value
     INTEGER(kind=4):: ll
     LOGICAL :: found

     ! begin
     IF( PRESENT(posin) )THEN
       iword = posin
       found = .TRUE.
     ELSE
       found = exists(kword,iword)   !see if key-word exists
       IF( PRESENT(posout) ) posout = iword
     END IF

     IF(found) THEN
       ival = INT(param(iword))    !associated param
       SELECT CASE (ival)
       CASE (:-1)                    !full label
         get_name = names(-ival)
       CASE (0)                      !label is the next word
         IF(nwopa > iword )THEN
           get_name = TRIM(words(iword+1)(1:MIN(mstl,mnam)))//'                        '
         ELSE
           get_name = TRIM(words(iword)(1:MIN(mstl,mnam)))//'                        '
         END IF
       CASE (1:)                     !label is a two digit number
         get_name = ''
         IF (PRESENT(stype)) THEN
           ll = MIN0(mnam-2,LEN_TRIM(stype))
           get_name(1:ll) = stype(1:ll)
         END IF
         get_name = TRIM(get_name)//TRIM(inttoch(ival,2))
       END SELECT

     ELSE IF ( PRESENT(texts) )THEN
       IF( texts(1:1) == '!')  THEN
         IF ( PRESENT(kword) )THEN
           WRITE(lures,1,ERR=9999) TRIM(kword),texts(2:LEN_TRIM(texts))
           WRITE(*    ,1,ERR=9999) TRIM(kword),texts(2:LEN_TRIM(texts))
         ELSE
           WRITE(lures,1,ERR=9999) TRIM(kword),texts(2:LEN_TRIM(texts))
           WRITE(*    ,1,ERR=9999) TRIM(kword),texts(2:LEN_TRIM(texts))
         END IF
         CALL runend('GET_NAME: COMPULSORY PARAM. NOT FOUND')
       ELSE IF ( PRESENT(defau) )THEN
         get_name = defau
       END IF

     END IF
     IF( PRESENT(texts) ) THEN
       IF( texts(2:2) /= '-' )THEN
         IF ( PRESENT(kword) )THEN
           WRITE (lures,"(9X,A,'.... ',A,' = ',(A))",ERR=9999) texts(2:LEN_TRIM(texts)),    &
             TRIM(kword),TRIM(get_name)
         ELSE
           WRITE (lures,"(9X,A,'........... = ',(A))",ERR=9999) texts(2:LEN_TRIM(texts)),   &
             TRIM(get_name)
         END IF
       END IF
     END IF
     IF(PRESENT(founds)) founds = found

   RETURN
   ! format for error messagge
  1 FORMAT(/,4X,'*** ERROR: ',A,' = ',A,/,15X,'IS A COMPULSORY PARAMETER. SPECIFY IT !')
  9999 CALL runen2('')
   END FUNCTION get_name

   INTEGER(kind=4) FUNCTION lfile(name)
   !-----------------------------------------------
   !     Returns the original leng of file name
   !-----------------------------------------------
   IMPLICIT NONE

     !Dummy arguments
     CHARACTER(len=*),INTENT(IN):: name
     ! local variables

     lfile = INDEX(name,'.@',BACK=.TRUE.) - 1
     IF (lfile < 0) lfile=LEN_TRIM(name)

   RETURN
   END FUNCTION lfile

   FUNCTION copy_file( old_name, new_name ) RESULT(dummy)
!#ifdef __INTEL_COMPILER
  USE ifport
!#else
!  INTEGER :: SYSTEM
!#endif
     CHARACTER(*), INTENT(IN) :: old_name, new_name
     CHARACTER(400) :: command
     INTEGER :: dummy

!#if UNIX
!     command = 'cp -f '//old_name//' '//new_name//' '
!     dummy = SYSTEM(TRIM(command))
!#else
!#  ifdef __INTEL_COMPILER
     command = 'COPY '//old_name//' '//new_name
     dummy = SYSTEM(TRIM(command))
!#  else
!     command = 'COPY '//old_name//' '//new_name
!     CALL SYSTEM(TRIM(command))
!     dummy = 0
!#  endif
!#endif
  END FUNCTION copy_file

  FUNCTION if_exist_copy_file( old_name, new_name ) RESULT(dummy)
!#ifdef __INTEL_COMPILER
  USE ifport
!#else
!  INTEGER :: SYSTEM
!#endif
     CHARACTER(*), INTENT(IN) :: old_name, new_name
     CHARACTER(400) :: command
     INTEGER :: dummy

!#if UNIX
!     command = 'if ( test -e '//old_name//' ); then '//    &
!               ' cp -f '//old_name//' '//new_name//' ; fi'
!     dummy = SYSTEM(TRIM(command))
!#else
!#  ifdef __INTEL_COMPILER
     command = 'IF EXIST '//old_name//' '//    &
               'COPY '//old_name//' '//new_name
     dummy = SYSTEM(TRIM(command))
!#  else
!     command = 'IF EXIST '//old_name//' '//    &
!               'COPY '//old_name//' '//new_name
!     call SYSTEM(TRIM(command))
!     dummy = 0
!#  endif
!#endif
  END FUNCTION if_exist_copy_file


  FUNCTION if_exist_remove_file( file_name1, file_name2 ) RESULT(dummy)
!#ifdef __INTEL_COMPILER
  USE ifport
!#else
!#if UNIX
!  INTEGER :: SYSTEM
!#else
!  USE dflib
!#endif
!#endif
     CHARACTER(*), INTENT(IN) :: file_name1, file_name2
     CHARACTER(400) :: command
     INTEGER :: dummy

!#if UNIX
!     command = 'if ( test -e '//file_name1//' ); then '//    &
!               ' rm -f '//file_name2//' ; fi'
!     dummy = SYSTEM(TRIM(command))
!#else
!#ifdef __INTEL_COMPILER
     command = 'IF EXIST '//file_name1//' DEL '//file_name2
     dummy = SYSTEM(TRIM(command))
!#else
!     command = 'IF EXIST '//file_name1//' DEL '//file_name2
!     CALL SYSTEM(TRIM(command))
!     dummy = 0
!#endif
!#endif
  END FUNCTION if_exist_remove_file


    SUBROUTINE rdtitle(ptype,istra,text)
    IMPLICIT NONE

    INTEGER (kind=4),INTENT(IN) :: istra
    CHARACTER(len=midn),INTENT(IN) ::  ptype
    CHARACTER(len=mttl),INTENT(OUT) ::  text

    CHARACTER(len=mich):: inttoch

      CALL listen('RTITLE')   !read TITLE

      IF(words(1) == 'TITLE' ) THEN
        text = ''
        WRITE(lures,"(//,10x,a60,//)",err=9999) card(1:60)
        text = TRIM(card)
      ELSE  !add a default title
        backs = .TRUE.
        IF (ptype == 'IMPLSP') THEN
          text = 'IMPLICIT SPRINGBACK (STRATEGY '//TRIM(inttoch(istra,2))//')'
        ELSE IF (ptype == 'CUTTIN') THEN
          text = 'CUTTING (STRATEGY '//TRIM(inttoch(istra,2))//')'
        ELSE IF (ptype == 'TRAROT') THEN
          text = 'POSITIONING (STRATEGY '//TRIM(inttoch(istra,2))//')'
        END IF
      END IF

    RETURN
    9999 CALL runen2('')
    END SUBROUTINE rdtitle
 END MODULE c_input
