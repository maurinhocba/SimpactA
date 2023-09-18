      SUBROUTINE scale_mass( mscal )
      !********************************************************************
      !
      !***   evaluation mass scale factor
      !
      !********************************************************************
      USE ctrl_db, ONLY : dtscal,dtuser,tdtuse,tdtsca,therm
      IMPLICIT NONE
       
      REAL(kind=8), INTENT(OUT) :: mscal  
 
      REAL (kind=8) :: deltc,tdelt,mtime,thtim
      INTERFACE
        INCLUDE 'elemnt.h'
      END INTERFACE

      !***  loop over all the sets

      deltc = 1.d0        !set a initial value
      CALL elemnt ('INCDLT', deltc= deltc) !compute elements Maximum Value
      IF(dtuser == 0d0)THEN        !check if the user gave a Maximum value
        mtime = deltc*dtscal       !use Maximum value
      ELSE
        mtime = MIN(deltc*dtscal,dtuser) !compare with user value
      END IF
          
      tdelt = 1.d0
      CALL elemnt ('TINCDT', deltc= tdelt)
      IF(tdtuse == 0d0)THEN
        thtim = tdelt*tdtsca
      ELSE
        thtim = MIN(tdelt*tdtsca,tdtuse)
      END IF
      
      mscal = thtim/mtime ! mass scale factor 

      RETURN
      END SUBROUTINE scale_mass
