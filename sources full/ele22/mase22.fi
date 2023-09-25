 SUBROUTINE mase22(nnode,nelem,matno,lnods,sname)

 IMPLICIT NONE
 INTEGER (kind=4) nnode,nelem,matno(nelem),lnods(nnode,nelem)
 CHARACTER (len=*) :: sname

 INTEGER (kind=4) ielem,inode

 WRITE(17,ERR=9999) 10,nelem,sname
 WRITE(17,ERR=9999) nnode,1,0
 WRITE(17,ERR=9999) matno(1),((lnods(inode,ielem),inode=1,nnode),ielem=1,nelem)

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE mase22
