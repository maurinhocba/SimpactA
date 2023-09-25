SUBROUTINE mshcri (istep, ttime, dtime, actmsh, endmsh)
!-----------------------------------------------------------------------
!  Verify mesh modification criteria
!-----------------------------------------------------------------------
USE meshmo_db
IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4),INTENT(IN):: istep !present step                               
  REAL(kind=8),INTENT(IN):: ttime, & !present time
                            dtime    !present time increment
  LOGICAL,INTENT(OUT):: actmsh, & !.TRUE. if mesh modifications required in step
                        endmsh    !.TRUE. if mesh modifications requires at end step 

  !--- Local variables  
  REAL (kind=8) :: toler,time1
  INTEGER (kind=4) :: i, is
  LOGICAL :: found

  INTERFACE
    INCLUDE 'elemnt.h'
  END INTERFACE

  toler = dtime*(1.0000000001d0)/2d0        !tolerance to check
  actmsh = .FALSE. !initializes
  endmsh = .FALSE. !
  
  IF (r_nrfset > 0) THEN    !for remeshing criteria
  
    DO is=1,r_nrfset
      r_meshmo(is) = .FALSE.   !initializes
      r_z_actv     = .FALSE.   !initializes
      IF( .NOT.r_last )THEN
        !initializes using step frequency criterion
        r_last =  MOD(istep-r_s_star, r_s_freq) == 0
        !check time frequency criterion
        time1 = MODULO(ttime-r_t_star,r_t_freq)
        IF(time1 < toler .OR. time1+toler > r_t_freq ) r_last = .TRUE.
        !check specific time criterion
        DO i=1,r_size
          IF( ABS(ttime-r_times(i)) < toler ) r_last = .TRUE.
        END DO
        !now check the second part
        IF( r_last )THEN  !
          !check criteria for remeshing
          IF(r_crit == 0)THEN
            !IF( r_min_dt_r * ini_dtime > dtime ) meshmo = 1
            !IF( r_min_size > p_elm_size ) meshmo = 1
            r_meshmo(is) = .TRUE.
          ELSE
            !diferent element distortion criteria for TR2D only
            CALL elemnt('ELMDIS',name=r_refset(is),flag2=found)
            r_meshmo(is) = r_z_actv  !.TRUE. -> if elements are distorted
          END IF
        END IF
      ELSE
        r_last = .FALSE.
      END IF
    END DO
  END IF

  !actmsh = .TRUE. for only mesh modification in step 
  IF( ASSOCIATED(r_meshmo) .AND. .NOT.endmsh ) actmsh = ANY(r_meshmo)

RETURN
END SUBROUTINE mshcri
