 SUBROUTINE actcd( )

 !keeps the coordinates at the begining of the stage ==> COORS

 USE npo_db, ONLY : coora,coors,coorc,psia,psic
 USE ctrl_db, ONLY: ndime, npoin,ndofn, initial_displacements
 IMPLICIT NONE

 !Local variable
 INTEGER (kind=4) :: i

 IF( initial_displacements )THEN
   DO i=1,npoin
     coora(:,i) = coorc(:,i) + coora(:,i)  !incremental displacement at the begining of the strategy
   END DO
   IF( ndofn == 8 )THEN
     DO i=1,npoin
       psic(:,i) = psic(:,i) + psia(:,i)   !incremental displacement at the begining of the strategy
     END DO
   END IF
   initial_displacements = .FALSE.
 END IF
 IF (ASSOCIATED(coors)) DEALLOCATE(coors)
 ALLOCATE(coors(ndime,npoin))
 DO i=1,npoin
   coors(:,i) = coora(:,i)
   coorc(:,i) = coora(:,i)  !no incremental displacement at the begining of the strategy
   IF( ndofn == 8 ) psia(:,i) = psic(:,i)      !add displacement at the begining of the strategy
 END DO
 RETURN
 END SUBROUTINE actcd
