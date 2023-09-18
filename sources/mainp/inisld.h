 SUBROUTINE inisld(npsdf,nesdf,ftsdf,j,np)
 !**********************************************************************
 !
 !     APPLIES and generates DATA for sliding nodes
 !     for 3-D problems
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: np    !size (to check it is not exceeded)
 INTEGER (kind=4),INTENT(IN OUT) :: j,npsdf(np),nesdf(4*np)
 REAL (kind=8), INTENT(IN OUT) :: ftsdf(4*np)
 END SUBROUTINE inisld
