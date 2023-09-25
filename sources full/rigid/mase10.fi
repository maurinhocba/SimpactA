 SUBROUTINE mase10(nnode,nelem,matno,lnods,sname,ntype,nmast)

 IMPLICIT NONE
 INTEGER (kind=4), POINTER :: matno(:)
 INTEGER (kind=4), INTENT(IN) :: nnode,nelem, &
                   lnods(nnode,nelem),ntype,nmast
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) ielem,inode,matn

 IF( nelem == 0 .OR. nnode == 0 ) RETURN
 IF( ASSOCIATED(matno) )THEN
   matn = matno(1)
 ELSE
   matn = 0
 END IF
 WRITE(17,ERR=9999) 10,nelem,sname
 WRITE(17,ERR=9999)  nnode,ntype,nmast
 WRITE(17,ERR=9999) matn,((lnods(inode,ielem),inode=1,nnode),ielem=1,nelem)


 RETURN
 9999 CALL runen2('')
 END SUBROUTINE mase10
