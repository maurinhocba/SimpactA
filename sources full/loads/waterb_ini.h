 SUBROUTINE waterb_ini (numfl,ndime,flparm,headf)
 !Reads follower loads
 USE loa_db,ONLY: foll_seg, flpar
 !USE npo_db,ONLY: coora
 IMPLICIT NONE

   !--- Dummy variables
   INTEGER(kind=4):: numfl, ndime
   TYPE(flpar):: flparm
   TYPE(foll_seg),POINTER:: headf
 END SUBROUTINE waterb_ini
