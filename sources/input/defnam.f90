 SUBROUTINE defnam(filnam,nstra)
 !
 !  set output file name-root according to problem
 !
 USE meshmo_db, ONLY: nrmstra !type of mesh modification active
 USE param_db, ONLY: mich
 USE c_input, ONLY : lfile
 IMPLICIT NONE

   !Dummy variables
   CHARACTER(len=*),INTENT(IN OUT):: filnam  !initial output file name root
   INTEGER(kind=4):: nstra ! present strategy
   !Functions
   CHARACTER(len=mich):: inttoch
   !Local variables
   INTEGER(kind=4):: lno

   lno = lfile(filnam) !original length ofr output file name-root
   filnam = filnam(1:lno)//'.@'//TRIM(inttoch(nstra,0))  !add strategy
   !  add suffix for mesh modifications
   IF (nrmstra > 0) filnam=TRIM(filnam)//'_rm'//TRIM(inttoch(nrmstra,0))

 RETURN
 END SUBROUTINE defnam
