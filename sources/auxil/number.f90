INTEGER(kind=4) FUNCTION number(string)
IMPLICIT NONE

  !--- Dummy variables
  CHARACTER(len=*),INTENT(IN):: string
  !--- Local variables
  INTEGER(kind=4):: i, io, n, nchar, signo
  CHARACTER(len=1):: ch

  number = 0
  nchar = LEN_TRIM(string)
  IF (nchar == 0) CALL runen2('Empty number string.')

  io = 0
  signo = 1
  DO i=1,nchar
    io = io + 1
    IF (string(i:i) == ' ') CYCLE
    IF (string(io:io) == '-') THEN
      signo = -1
      io = io + 1
    ELSE IF (string(io:io) == '+') THEN
      io = io + 1
    END IF
    EXIT
  END DO

  DO i=io,nchar
    ch = string(i:i)
    IF (ch == ' ') EXIT
    n = ICHAR(ch)-48
    IF (n >= 0 .AND. n <= 9) THEN
      number = number*10 + n
    ELSE
      CALL runen2('String cannot be converted into an integer.')
    END IF
  END DO
  number = signo*number

RETURN
END FUNCTION number
