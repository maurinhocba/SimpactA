  SUBROUTINE runend(message)
  !*************************************************************************
  !
  !****this routine stops the run and prints the last read card
  !    the message MESSAGE is output to different places (screen, debug, report)
  !
  !*************************************************************************
  USE lispa0
  USE name_db, ONLY: gen_data_file
  IMPLICIT NONE

    !--- Dummy variables
    CHARACTER(len=*),INTENT(IN):: message
    !--- Local variables
    INTEGER(kind=4):: i

    WRITE(55,"(//,5X,34('* '),//,25X,'AN ERROR HAS BEEN DETECTED:',/)",ERR=9999)
    WRITE(55,*,ERR=9999) TRIM(message)
    WRITE(55,"(5X,'LAST CARD READ : ',/,A,//,5X,34('* '))",ERR=9999) TRIM(card)
    WRITE(55,"(5X,'Number of words: ',i2,'  Number of pars :',i2)",ERR=9999) nnwor,nnpar
    WRITE(55,"(5X,A,' = ',e20.10)",ERR=9999) (TRIM(words(i)),param(i),i=1,nwopa)

    WRITE(*,"('* * * AN ERROR HAS BEEN DETECTED * * *')")
    IF( gen_data_file )THEN
      CLOSE(ludat,STATUS='DELETE')   !delete input data file
    ELSE
      CLOSE(ludat,STATUS='KEEP')   !delete input data file
    END IF

    STOP '* * * See detailes in *.rep file * * *'
    9999 CALL runen2('')
  END SUBROUTINE runend

 SUBROUTINE runen2(message)
 !*************************************************************************
 !
 !****this routine stops the run when a WRITE error is detected
 !    the message MESSAGE is output to different places (screen, debug, report)
 !
 !*************************************************************************
 USE lispa0
 USE name_db, ONLY: gen_data_file
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*), INTENT(IN) :: message

   IF(LEN_TRIM(message) > 0) WRITE(*,1) TRIM(message)   !print message to screen
   IF( gen_data_file )THEN
     CLOSE(ludat,STATUS='DELETE')   !delete input data file
   ELSE
     CLOSE(ludat,STATUS='KEEP')   !delete input data file
   END IF
   WRITE(uf,1) TRIM(message)      !print message to report file
   WRITE(55,1) TRIM(message)      !print message to report file
   CLOSE (55)
   STOP                           !stop execution

 1 FORMAT(//,'  AN ERROR HAS BEEN DETECTED:',                          &
           /,'  FAILED WRITTING, THE DISK IS POSSIBLY FULL.',          &
          //,'  PLEASE, CHECK AND FREE DISK SPACE BEFORE CONTINUING.', &
          //,A)
 END SUBROUTINE runen2

 SUBROUTINE runen3(message)
 !*************************************************************************
 !
 !****this routine stops the run and prints an error message
 !
 !*************************************************************************
 USE lispa0
 USE outp_db, ONLY : time
 USE name_db, ONLY: gen_data_file
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*), INTENT(IN) :: message
   INTERFACE
     INCLUDE 'byebye.h'
   END INTERFACE

   WRITE(6,"(//,'* * * AN ERROR HAS BEEN DETECTED * * *')")
   WRITE(6,"(A)") 'ERROR: '//TRIM(message)
   WRITE(55,"(//,'* * * AN ERROR HAS BEEN DETECTED * * *')",ERR=9999)
   WRITE(55,*,ERR=9999) 'ERROR: '//TRIM(message)
   IF( gen_data_file )THEN
     CLOSE(ludat,STATUS='DELETE')   !delete input data file
   ELSE
     CLOSE(ludat,STATUS='KEEP')   !delete input data file
   END IF
   CALL timing(19,2)
   CALL byebye(time)
   CLOSE(lures,STATUS='KEEP') !
   CLOSE (55)

   STOP
   9999 CALL runen2('')
 END SUBROUTINE runen3
