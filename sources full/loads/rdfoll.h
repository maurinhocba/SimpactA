 SUBROUTINE rdfoll (numfl,ndime,iwrit,fltype,flparm,headf,tailf,hydfl,curv_name, &
                    fluid)
 !Reads follower loads
 !USE c_input,ONLY: ludat, lures, listen, exists, getint, getrea, backs
 USE loa_db !,ONLY: foll_seg, flpar, hydf_load, ini_hydfl, rdhyd, new_foll, add_foll
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*):: fltype   !len=6
   CHARACTER(len=*):: curv_name
   INTEGER(kind=4):: numfl, ndime, iwrit
   TYPE(flpar):: flparm
   TYPE(foll_seg),POINTER:: headf, tailf
   TYPE(hydf_load),POINTER:: hydfl
   LOGICAL :: fluid

 END SUBROUTINE rdfoll
