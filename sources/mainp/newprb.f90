  SUBROUTINE newprb()
 !********************************************************************
 !
 !*** set output file name root
 !
 !********************************************************************
 USE name_db,ONLY: output, rsfprj, rsfstg
 USE c_input, ONLY : lfile,rdtitle
 USE ctrl_db, ONLY : nstra,text,ptype
 IMPLICIT NONE

   !---Local variables
   INTEGER(kind=4):: lrn     !lenght of the original root name (variable)

   CALL defnam(output,nstra)   !determines the output file name root according to mesh modifications
   lrn = lfile(output)   !length of original root name
   rsfprj = output(1:lrn)                    !name of the project (original output name)
   rsfstg = output(lrn+1:LEN_TRIM(output))   !strategy suffix

   CALL rdtitle(ptype,nstra,text)

 RETURN
 END SUBROUTINE newprb
