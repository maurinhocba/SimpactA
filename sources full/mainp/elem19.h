 SUBROUTINE elem19(TASK, elsnam, dtime, ttime, istop, lnod, flag1, flag2)

 !master routine for element 19 (TLF) 2-D solid triangular element with thermal coupling

 !USE ctrl_db, ONLY: ndofn, npoin, ntype
 !USE outp_db, ONLY: sumat, iwrit
 USE ele19_db
 !USE npo_db
 !USE sms_db, ONLY : selective_mass_scaling, sms_ns , sms_name, sms_thl, sms_alp
 IMPLICIT NONE

 !--- Dummy variables
 CHARACTER(len=*),INTENT(IN):: TASK
 ! optional parameters
 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL :: elsnam
 INTEGER (kind=4), OPTIONAL :: istop
 INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
 REAL (kind=8), OPTIONAL :: dtime,ttime

 END SUBROUTINE elem19
