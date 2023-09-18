 SUBROUTINE ubicmx ( ndofn,nelem,maxnn,neq,maxa,ifpre,lnods,  &
                     maxav,npsdf,nesdf)
 !**********************************************************************
 !
 !     routine to find the position of the diagonal elements in matrix a
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: ndofn,maxnn,nelem,neq,        &
                  ifpre(:,:),npsdf(:),nesdf(:),lnods(:,:)
 INTEGER (kind=4),INTENT(OUT) :: maxa,maxav(:)
 END SUBROUTINE ubicmx
