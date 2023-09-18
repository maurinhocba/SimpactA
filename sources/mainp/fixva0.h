       SUBROUTINE fixva0 (iwrit,ndime,ndofn,actio,rot_free,nvelr)

  !***  READ the fixed values.

       IMPLICIT NONE
       LOGICAL :: rot_free
       CHARACTER(len=*),INTENT(IN):: actio
       INTEGER(kind=4) :: iwrit,ndime,ndofn,nvelr

       END SUBROUTINE fixva0
