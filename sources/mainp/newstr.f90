 SUBROUTINE newstr(actio)
 !********************************************************************
 !
 !*** CHECK BEGINNING OF A NEW STRATEGY
 !    ACTIO = NSTRA0 or MESHMO
 !
 !********************************************************************
 USE c_input,ONLY: listen, exists, del_old_files, lfile, rdtitle
 USE name_db,ONLY: output, rsfstg
 USE ctrl_db, ONLY : nstra,text,ptype
 USE meshmo_db,ONLY: nrmstra, ntrmstra, r_meshmo
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*),INTENT(IN OUT):: actio  !(IN) NSTRA0 or MESHMO or ENDMSH
   !--- Local variables
   INTEGER(kind=4):: lrn

   IF (TRIM(actio) /= 'MESHMO') THEN  !actio = 'NSTRA0' or 'ENDMSH'
     CALL listen('NEWSTR')            !read a line
     IF (.NOT.exists('NEWSTR').OR.exists('ENDOFD')) actio='STOP'  !error or end of data
   END IF

   IF (TRIM(actio) /= 'STOP') THEN    !NEW_STR card read OK
     IF (exists('IMPLSP')) THEN       !IMPLicit SPringback
       ptype = 'IMPLSP'
     ELSE IF (exists('POSITI')) THEN  !POSITIoning of a tool or mesh
       WRITE(6,"(//,13x,'*** MOVING NODES AND INITIALIZING NEW STAGE ANALYSIS ***',/)")
       ptype = 'TRAROT'               !TRAnslate and/or ROTate
     ELSE IF (TRIM(actio) == 'MESHMO') THEN  !MESH MOdifications
       WRITE(6,"(/,13x,'*** MODIFYING MESH AND INITIALIZING NEW STAGE ANALYSIS ***',/)")
     ELSE IF (TRIM(actio) == 'ENDMSH') THEN  !MESH MOdifications
       WRITE(6,"(/,13x,'*** MODIFYING MESH AT END STAGE ANALYSIS ***',/)")
     ELSE                             !default strategy
       ptype = 'EXPLCT'               !EXPLiCit time integration
       WRITE(6,"(//,19x,'*** INITIALIZING NEW STAGE ANALYSIS ***',/)")
     END IF
     !     For MESH MOdifications in STEP
     IF (TRIM(actio) == 'MESHMO') THEN
       IF (ANY(r_meshmo)) THEN   !Remesh is working
         nrmstra = nrmstra + 1     !Increase strategy remesh counter
         ntrmstra = ntrmstra + 1   !Increase total remesh counter
       END IF
     ELSE ! only if strategy has data
       nstra = nstra + 1       !number of the new strategy
       !  initializes number of mesh modifications for the strategy
       nrmstra = 0
       CALL rdtitle(ptype,nstra,text)
     END IF

     CALL defnam(output,nstra)   !Determines the output file name root
     lrn = lfile(output)   !length of the original output file name root
     rsfstg = output(lrn+1:LEN_TRIM(output))  !Name of the strategy (eg. output file: <rsfprj><rsfstg>.f16)
     CALL del_old_files(2) !erase files from previous run (if files exists)

   END IF

 RETURN
 END SUBROUTINE newstr
