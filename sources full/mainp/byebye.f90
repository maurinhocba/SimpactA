      SUBROUTINE byebye (time)
!*********************************************************************
!     writes cpu times
!*********************************************************************
      USE c_input, ONLY : openfi
      USE gvar_db, ONLY : actchk
      IMPLICIT NONE
      REAL (kind=8),INTENT(IN OUT) :: time(:)

      INTEGER (kind=4) :: i,iha,imi,isg
      REAL (kind=8) :: timeo,timeop

      IF(actchk)RETURN

      DO i=2,10
        time(i+20) = time(i)*100d0/time(1)
      END DO

      DO i=11,17
        time(i+20) = time(i)*100d0/time(19)
      END DO
      time(21) = time( 1)*100d0/time(19)
      time(32) = time(12)*100d0/time(19)
      time(38) = time(18)*100d0/time(19)
      timeo = time(19) - time(1) - time(12) - time(18)
      timeop = timeo*100d0/time(19)
      CALL openfi(2)
      WRITE (2,40,ERR=9999) time( 1),time(21),time( 2),time(22),&
     &                      time( 5),time(25),time( 6),time(26),&
     &                      time( 4),time(24),&
     &                      time( 7),time(27),time( 8),time(28),&
     &                      time( 9),time(29),time(10),time(30),&
     &                      time(12),time(32),   timeo,  timeop,&
     &                      time(13),time(33),time(14),time(34),&
     &                      time(15),time(35),time(16),time(36),&
     &                      time(17),time(37),time(11),time(31),&
     &                      time(18),time(38)
      WRITE (2,41,ERR=9999) time(19)
      iha = INT(             time(19)/3600.0              )
      imi = INT(    (time(19)-REAL(iha)*3600.0)/60.0      )
      isg = INT( time(19)-REAL(iha)*3600.0-REAL(imi)*60.0 )
      WRITE (2,42,ERR=9999) iha,imi,isg
      CLOSE(2,STATUS='keep')
      RETURN
   40 FORMAT(///,'T I M I N G   I N F O R M A T I O N',//,&
     & 5X,'ACCION                           CPU(SEC)         %',//,&
     & 5X,'INITIAL COMPUTATIONS .........',E14.4,2X,F6.2,' %',//,&     ! 1
     &15X,'RESTART INPUT-OUTPUT .........',E14.4,2X,F6.2,' %',//,&     ! 2
     &15X,'GAUSSIN VARIABLES INTIALIZAT..',E14.4,2X,F6.2,' %',//,&     ! 5
     &15X,'DATA INPUT ...................',E14.4,2X,F6.2,' %',//,&     ! 6
     &15X,'CONTACT DATA INPUT ...........',E14.4,2X,F6.2,' %',//,&     ! 4
     &15X,'INITIAL CONDITIONS ...........',E14.4,2X,F6.2,' %',//,&     ! 7
     &15X,'INITIAL STRESS STATE .........',E14.4,2X,F6.2,' %',//,&     ! 8
     &15X,'NODAL LOAD VECTOR COMPUT. ....',E14.4,2X,F6.2,' %',//,&     ! 9
     &15X,'LUMPED MASS COMPUTATION ......',E14.4,2X,F6.2,' %',//,&     !10
     &15X,'DUMPING OF SPRINGBACK DATA....',E14.4,2X,F6.2,' %',//,&     !12
     & 5X,'DYNAMIC ANALYSIS .............',E14.4,2X,F6.2,' %',//,&     !timeo
     &15X,'TIME INTEGRATION .............',E14.4,2X,F6.2,' %',//,&     !13
     &15X,'RESIDUAL FORCES COMPUTATION ..',E14.4,2X,F6.2,' %',//,&     !14
     &15X,'RESULTS OUTPUT ...............',E14.4,2X,F6.2,' %',//,&     !15
     &15X,'CONTACT COMPUTATIONS .........',E14.4,2X,F6.2,' %',//,&     !16
     &15X,'TIME INCREMENT EVALUATION ....',E14.4,2X,F6.2,' %',//,&     !17
     &15X,'MESH REFINAMENT EVALUATION ...',E14.4,2X,F6.2,' %',//,&     !11
     & 5X,'FILE CLOSING .................',E14.4,2X,F6.2,' %',//)     !18
   41 FORMAT(5X,'T O T A L                    ',E14.4,//)
   42 FORMAT(////,15X,'HOURS: ',I6,8X,'MINUTES:',I2,8X,'SECONDS:',I2,//)
 9999 CALL runen2('')
      END SUBROUTINE byebye
