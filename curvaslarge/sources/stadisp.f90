      SUBROUTINE stadisp(udf,ulf,nvrbl,nreq,nprq,empty,option,varname)
!=======================================================================
!     print Info data-input for different variables
!=======================================================================
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: udf,      & !file to read information (binary)
                                     ulf,      & !file write information (ASCII)
                                     nvrbl,    & !number of variables per node
                                     nreq,     & !number of nodes in the list
                                     nprq(nreq)  !list of nodes
      LOGICAL,INTENT(INOUT) :: empty             !.TRUE. if file ULF is empty
      CHARACTER (len=1) :: option                !character option
      CHARACTER (len=*) :: varname               !name of the variable

      !Local variables
      INTEGER (kind=4) :: j, ist, nstep, mm
      REAL (kind=8) :: aa
      !CHARACTER (len=1) :: postype

      IF (nreq == 0 .OR. nvrbl == 0) RETURN      !no nodes in the list

      !READ(udf,IOSTAT=ist) postype               !this is not a variable
      !IF (ist/=0) STOP ' File contain incorrect information.'

      nstep = 0
      DO  !loop to read all the steps
        READ(udf,IOSTAT=ist) aa      !read one step information
        IF (ist /= 0) EXIT           !if error EXIT
        nstep = nstep + 1
      END DO
      IF (nstep == 0) RETURN         !no step read Exit routine

      empty = .FALSE.                !modify flag

      ! writes information into file ULF=3
      WRITE(ulf,"(A,A)",err=9999) 'BEGIN_',TRIM(varname)  !first heading
      WRITE(ulf,"(A,A)",err=9999) 'OPTION= ',option       !second heading
      DO j=1,nreq,10
        mm = MIN(nreq,j+9)
        WRITE(ulf,"(A,10(I6,1X))",err=9999) 'POINTS= ',nprq(j:mm)
      END DO
      WRITE(ulf,"(A,I3)",err=9999) 'NCOMP= ',nvrbl
      WRITE(ulf,"(A,A,/)",err=9999) 'END_',TRIM(varname)

      RETURN
 9999 WRITE(6,"(//,A,/,A)")' AN ERROR HAS BEEN DETECTED WHILE'//         &
     &  ' WRITING TO DISK. THE DISK IS POSSIBLY FULL.',' PLEASE'//       &
     &  ' CHECK AND FREE DISK SPACE BEFORE CONTINUING.'
      STOP


      END SUBROUTINE stadisp
