 SUBROUTINE loadfl (ttime)

 ! calculates and applies follower load

 USE ctrl_db, ONLY: ndime,nload,dtime,ntype,npoin,endtm
 USE loa_db
 USE npo_db, ONLY : coora,resid,loass
 USE hydf_db
 IMPLICIT NONE
 REAL (kind=8), INTENT(IN) :: ttime
 INTERFACE
    INCLUDE 'aplfol.h'
 END INTERFACE

 ! Local
 INTEGER (kind=4) :: il
 TYPE (loa_set), POINTER :: loas
 !INTEGER, SAVE :: itera = 0  !to write debug values for waterbags
 !---------------------------------------------------------------


 loas => headls
 DO il = 1,nload ! loop over load sets
   ! if set contains follower load data
   IF (loas%numfl > 0) CALL aplfol(ndime,coora,resid,loass(il), &
                      dtime,ttime,ntype,loas%fltype,       &
                      loas%flpar,loas%headf, loas%factor)
   !IF(MOD(itera,100) == 0 .AND. TRIM(loas%fltype) == 'WATERB') &
   !  WRITE(55,"(3e15.5)")loas%flpar%press,loas%flpar%pext,loas%flpar%vol
   IF (ASSOCIATED (loas%hydfl) ) & ! if set contains hydroforming data
     CALL hydflo (loas%hydfl, loass(il), loas%factor)
   loas => loas%next
 END DO
 !itera = itera + 1

 END SUBROUTINE loadfl
