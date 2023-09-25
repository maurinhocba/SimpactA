 SUBROUTINE dump10(sname, nelem, nnode, rbnod, nmast, ntype, &
                   lnods, matno, tmass,                      &
                   heat, ngaus, shape, cartd, dvolu )
 !******************************************************************
 !
 !*** dumps an element set for restart
 !
 !******************************************************************
 IMPLICIT NONE

 CHARACTER (len=*) :: sname
 INTEGER (kind=4) :: nelem, nnode, rbnod, nmast, ntype,  &
                     lnods(:,:), ngaus
 INTEGER (kind=4), POINTER :: matno(:)
 REAL(kind=8) :: tmass
 LOGICAL :: heat
 REAL(kind=8), POINTER :: shape(:,:),cartd(:,:,:,:),dvolu(:,:)

 INTEGER (kind=4) :: ielem,i,j
 CHARACTER(len=midn):: straux
 LOGICAL (kind=4) :: flag

 straux = 'ELSET'
 WRITE (50,ERR=9999) straux
 straux = 'RIGID'
 WRITE (50,ERR=9999) straux,sname
 flag = ASSOCIATED( matno )
 WRITE (50,ERR=9999) nelem
 WRITE (50,ERR=9999) nnode, rbnod, nmast, ntype, flag, tmass, heat, ngaus

 IF( flag ) THEN
   DO ielem=1,nelem
     WRITE (50,ERR=9999) lnods(1:nnode,ielem) , matno(ielem)
   END DO
 ELSE
   DO ielem=1,nelem
     WRITE (50,ERR=9999) lnods(1:nnode,ielem)
   END DO
 END IF
 IF( heat )THEN
   WRITE (50,ERR=9999) ((shape(i,j),i=1,nnode),j=1,ngaus)
   DO ielem=1,nelem
     WRITE (50,ERR=9999) cartd(:,:,:,ielem),dvolu(:,ielem)
   END DO
 END IF
 RETURN
 9999 CALL runen2('')
 END SUBROUTINE dump10
