SUBROUTINE rembln (string1)

  ! removing blanks from a string
  USE param_db,ONLY: mnam
  IMPLICIT NONE

  CHARACTER (len=*), INTENT(INOUT) :: string1

  ! local
  CHARACTER (len=mnam) :: string2
  INTEGER(kind=4) :: i, n
  
  !=============================

  string2 = string1
  string1 = ''
  n = 0
  DO i=1,LEN_TRIM(string2)
    IF (string2(i:i) ==' ') CYCLE
    n = n + 1
    string1(n:n) = string2(i:i) 
  ENDDO

  RETURN

END SUBROUTINE rembln
