      SUBROUTINE stacpvt(udf,ulf,iener,nreqe,empty,option,varname)
!=======================================================================
!  Info data-input for cynetic energy, potential energy, velocity  & dtime
!=======================================================================
      IMPLICIT NONE

      INTEGER (kind=4),INTENT(IN) :: udf,      & !file to read information (binary)
                                     ulf,      & !file write information (ASCII)
                                     iener,    & !1 if present
                                     nreqe       !number of values
      LOGICAL,INTENT(INOUT) :: empty             !.TRUE. if file ULF is empty
      CHARACTER (len=1) :: option                !character option
      CHARACTER (len=*) :: varname               !name of the variable

      !Local variables
      INTEGER (kind=4) :: ist, nstep
      REAL (kind=8) :: aa
      !CHARACTER (len=1) :: postype

      IF (iener == 0) RETURN

      !READ(udf,IOSTAT=ist) postype
      !IF (ist/=0) STOP ' File contain incorrect information.'

      nstep = 0
      DO
        READ(udf,IOSTAT=ist) aa
        IF (ist /= 0) EXIT
        nstep = nstep + 1
      END DO
      IF (nstep == 0) RETURN

      empty = .FALSE.

      WRITE(ulf,"(A,A)",err=9999) 'BEGIN_',TRIM(varname)  !first heading
      WRITE(ulf,"(A,A)",err=9999) 'OPTION= ',option       !second heading
      WRITE(ulf,"(A,I3)",err=9999) 'VARIAB= ',nreqe
      WRITE(ulf,"(A,A,/)",err=9999) 'END_',TRIM(varname)

      RETURN
 9999 WRITE(6,"(//,A,/,A)")' AN ERROR HAS BEEN DETECTED WHILE'//         &
     &  ' WRITING TO DISK. THE DISK IS POSSIBLY FULL.',' PLEASE'//       &
     &  ' CHECK AND FREE DISK SPACE BEFORE CONTINUING.'
      STOP

      END SUBROUTINE stacpvt
