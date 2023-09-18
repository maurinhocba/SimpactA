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
 !local variables
 INTEGER(kind=4) :: i,n
 REAL(kind=8) :: aux,tn(3,3)

     IF( ALL(.NOT.isidf))THEN     !no clamped sides
       DO i=1,nnode               !for each node
         n = lnods(i)             !node
         IF(bottom) coorb(:,n) = coorb(:,n) - t3*thnew !bottom surface
         IF(top   ) coort(:,n) = coort(:,n) + t3*thnew !top surface
         ifact(n) = ifact(n) + 2  !increase weigthing factor
       END DO
     ELSE                         !one or more clamped sides
       n = 0                      !initializes number of clamped sides
       DO i=1,nnode               !loop over each side
         IF( isidf(i) )THEN       !if side clamped
           n = n+1                !increase number of clamped sides
           tn(:,n) = phin(:,i)/2d0  !normal to plane of symmetry
         END IF
       END DO
       IF( n == 1)THEN            !for one clamped side (or symmetry plane)
         CALL vecuni(3,tn(1,1),aux)
         aux = DOT_PRODUCT(tn(:,1),t3)    !proyect normal over plane
         tn(:,3) = t3 - aux*tn(:,1)       !proyect normal over plane
       ELSE                       !for two clamped side (or symmetry planes)
         CALL vecpro(tn(1,1),tn(1,2),tn(1,3)) !normal to both planes
       END IF
       CALL vecuni(3,tn(1,3),aux)             !unit vector
       DO i=1,nnode               !for each node
         n = lnods(i)        !node
         IF(bottom) coorb(:,n) = coorb(:,n) - tn(:,3)*thnew !bottom surface
         IF(top   ) coort(:,n) = coort(:,n) + tn(:,3)*thnew !top surface
         ifact(n) = ifact(n) + 2  !increase weigthing factor
       END DO
     END IF
 RETURN
 END SUBROUTINE gentbs
