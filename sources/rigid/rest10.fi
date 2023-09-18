 SUBROUTINE rest10( nelem, nnode,lnods, matno, flag, ndime, &
                   heat, ngaus, shape, cartd, dvolu)
 !******************************************************************
 !
 !*** restores an element set at restart
 !
 !******************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: nelem,nnode,ndime,ngaus
 INTEGER (kind=4), POINTER :: lnods(:,:),matno(:)
 REAL (kind=8), POINTER :: shape(:,:),cartd(:,:,:,:),dvolu(:,:)
 LOGICAL (kind=4) :: flag,heat

 INTEGER (kind=4) :: ielem,i,j


 IF( flag ) THEN
   DO ielem=1,nelem
     READ (51) lnods(:,ielem) , matno(ielem)
   END DO
 ELSE
   DO ielem=1,nelem
     READ (51) lnods(:,ielem)
   END DO
 END IF
 IF( heat )THEN
   ALLOCATE( shape(nnode,ngaus),cartd(nnode,ndime,ngaus,nelem),dvolu(ngaus,nelem))
   READ (51) ((shape(i,j),i=1,nnode),j=1,ngaus)
   DO ielem=1,nelem
     READ (51) cartd(:,:,:,ielem),dvolu(:,ielem)
   END DO
 END IF

 RETURN
 END SUBROUTINE rest10
