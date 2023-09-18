 SUBROUTINE ubicmx ( ndofn,nelem,maxnn,neq,maxa,ifpre,lnods,  &
                     maxav,npsdf,nesdf)
 !**********************************************************************
 !
 !     routine to find the position of the diagonal elements in matrix a
 !
 !**********************************************************************
 USE lispa0
 USE kinc_db, ONLY : nn
 IMPLICIT NONE
 INTEGER (kind=4),INTENT(IN) :: ndofn,maxnn,nelem,neq,        &
                  ifpre(:,:),npsdf(:),nesdf(:),lnods(:,:)
 INTEGER (kind=4),INTENT(OUT) :: maxa,maxav(:)
 !
 INTEGER (kind=4) :: i,ii,j,k,l,n,mini,lb,le,jj

 !     initializes maxav
 DO i=1,neq+1
   maxav(i) = i
 END DO

 !    find the first element of each column

 DO n = 1,nelem
   mini = neq

 !       find the lowest equation number of each element  mini
   DO j = 1,maxnn
     i = lnods(j,n)
     IF(i > 0) THEN
       DO k = 1,ndofn
         ii = ifpre(k,i)
         IF(ii > 0) THEN
           IF(ii < mini) mini = ii
         ELSE IF(ii < 0 .AND. ii > -nn) THEN
           lb = npsdf(-ii)
           le = npsdf(-ii+1) - 1
           DO l = lb,le
             jj = nesdf(l)
             IF(jj > 0 .AND. jj < mini) mini  =  jj
           END DO
         END IF
       END DO
     END IF
   END DO

 !       compares mini with the previous values for each equation
   DO j = 1,maxnn
     i = lnods(j,n)
     IF(i > 0) THEN
       DO k = 1,ndofn
         ii = ifpre(k,i )
         IF(ii > 0 )THEN
           IF( maxav(ii) > mini) maxav(ii) = mini
         ELSE IF(ii < 0 .AND. ii > -nn) THEN
           lb = npsdf(-ii)
           le = npsdf(-ii+1) - 1
           DO l = lb,le
             jj = nesdf(l)
             IF(jj <= 0 )CYCLE
             IF( maxav(jj) > mini) maxav(jj) = mini
           END DO
         END IF
       END DO
     END IF
   END DO
 END DO

 !     determines the position of the diag. elements in vector stiff
 !MAXAV = 1 !this is to use a FULL SQUARE MATRIX

 maxa = 1
 DO i = 2,neq+1
   j = maxav(i)
   maxav(i) = maxav(i-1) + i - maxa
   maxa = j
 END DO
 maxa = maxav(neq+1)
 WRITE(lures,"(/,'  Number of active Equations        =  ',i10,/,       &
              &  '  Number of elements under Skyline  =  ',i10,/)")     &
                    neq,maxa
 RETURN

 END SUBROUTINE ubicmx
