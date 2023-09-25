 SUBROUTINE gentbs(isidf,lnods,bottom,top,coorb,coort,ifact,t3,thnew,phin,nnode)

 ! GENerates Top and Bottom Surfaces of the shell for contact analysis
 ! for three node shell elements (rotation free) LBST,CBST,NBST

 IMPLICIT NONE
 ! dummy arguments
 INTEGER(kind=4), INTENT(IN) :: nnode        !number of nodes
 INTEGER(kind=4), INTENT(IN) :: lnods(nnode) !connectivities
 LOGICAL, INTENT(IN)  :: isidf(nnode),     & !side boundary condition
                         bottom,top          !flags to indicate generation
 REAL(kind=8), INTENT(IN) :: t3(3),        & !element normal
                             thnew,        & !present element thickness
                             phin(3,nnode)   !symmetry plane normal
 INTEGER(kind=4), INTENT(IN OUT) :: ifact(:)
 REAL(kind=8), INTENT(IN OUT) :: coorb(:,:),coort(:,:)
 END SUBROUTINE gentbs
