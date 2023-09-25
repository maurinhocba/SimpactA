 SUBROUTINE scater(elvec,nrows,nnode,glvec,nglob,npoin,lnods)
 !***********************************************************************
 !
 !****this routine performs scatter operations
 !          elvec(nrows,nnode) ---> glvec(nrows,npoin)
 !
 !***********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) nrows,nnode,npoin,nglob,lnods(nnode)
 REAL (kind=8) elvec(nrows,nnode), glvec(nglob,npoin)

 INTEGER (kind=4)inode,lnode


 DO inode=1,nnode
   lnode=lnods(inode)
   glvec(1:nrows,lnode)=glvec(1:nrows,lnode)+elvec(1:nrows,inode)
 END DO

 END SUBROUTINE scater
