 SUBROUTINE import ( actio )
 !
 !   import element data set
 !
 USE param_db,ONLY: mnam
 USE c_input
 USE ctrl_db, ONLY: ndime, neulr, ndofn, npoin
 USE npo_db, ONLY : label,coord,coora,euler,ifpre,iffix,naeul
 USE mat_dba, ONLY : snn
 USE gvar_db, ONLY : fimpo,lab1,renum,seque,inter,overw
 USE meshmo_db, ONLY : nodset,numpo
 USE esets_db, ONLY: nelms, elsets, add_name, rot_free
 USE ifx_db, ONLY : nd1
 USE name_db, ONLY: output

 IMPLICIT NONE
 !dummy argument
 CHARACTER(len=*),INTENT(IN OUT):: actio !NEW, NSTRA0
 !Local variable
 REAL (kind=8) :: ang(3)
 CHARACTER(len=mnam) :: fname    !filename to read data
 LOGICAL :: found  !flag
 INTEGER (kind=4) :: nimpf  !Number of IMPort Files
 INTEGER (kind=4) :: i,j,nm,ns,itype,idime,ieulr,n
 INTEGER(kind=4), ALLOCATABLE :: lab0(:),   & !labels to renumber
                                 mn(:,:,:), & !materials labels (olds and news)
                                 sn(:,:,:)    !section labels (olds and news)
 LOGICAL, ALLOCATABLE :: impda(:,:)           !import flags
  !(1):keep/change name set             (2):keep/change node labels
  !(3):add fixed/sequential renumbering (4:8):cooro,matsc,bound,displ,inter
  !(9):add/overwrite

 CHARACTER(len=mnam), ALLOCATABLE :: els_name(:) !import set names

 INTERFACE
   INCLUDE 'elmdat.h'
   INCLUDE 'imcoor.h'
   INCLUDE 'immtsc.h'
   INCLUDE 'inrotm.h'
 END INTERFACE

 !-------
 CALL listen('IMPORT')
 IF( .NOT. exists('IMPORT' ))THEN  !if no import files
   backs = .TRUE.
   RETURN
 END IF
 nimpf=getint('IMPORT',20,' Number of imported sets ..........') !number of imported sets
 IF( nimpf == 0 )nimpf = 20  !maximum number of sets to import
 ! allocate auxiliar arrays
 ALLOCATE( lab0(nimpf), impda(9,nimpf), els_name(nimpf), mn(2,10,nimpf), sn(2,20,nimpf))
 i = 0                         !initializes counter of imported files
 DO
   CALL listen('IMPORT')       !read import file
   IF( exists('ENDIMP' ))THEN  !if list of import files ended
     IF( nimpf > i ) nimpf = i !actual number of imported sets
     EXIT                      !exit loop
   END IF
   i = i+1                   !update number of imported sets
   IF( i > nimpf )CALL runend('IMPORT: too many imported element sets')
   fname = get_name('FILE  ',found, '!Import File Name:')  !file name
   CALL openfi(70+i,fname)   !open file
   IF( exists('ELSNAM',j))THEN   !see if a new name is provided
     els_name(i) = get_name(posin=j,stype='ESET')  !file name
     impda(1,i) = .TRUE.  !remember that name have changed (unnecessary perhaps)
     READ(70+i)              !read a record
   ELSE
     READ(70+i)els_name(i)       !read the name from file
     impda(1,i) = .FALSE.        !original name
   END IF
   IF( exists('OVERWR',j) )THEN
     impda(9,i) = .TRUE.         !overwrite set with the same name
   ELSE
     impda(9,i) = .FALSE.        !do not overwrite set with the same name
   END IF
   ! if overwrite, name exists
   IF(.NOT.impda(9,1)) CALL add_name(els_name(i),2)       !add to element set names and check if exists
   impda(2,i) = exists('RENUMB') !check if node labels will be changed
   IF( exists('FROM  ',j) )THEN
     impda(3,i) = .TRUE.         !sequential renumbering
     lab0(i) = INT(param(j))-1   !first label (-1)
   ELSE IF( exists('ADD   ',j) )THEN
     impda(3,i) = .FALSE.        !same labels + a fixed value
     lab0(i) = INT(param(j))     !value to add
   ELSE
     CALL runend('IMPORT: ADD or FROM is compulsory')
   END IF
   READ(70+i)impda(4:8,i)     !read flags (cooro,matsc,bound,displ,inter) from binary file
 END DO
 !-------

 !------  READ coordinates and displacements of nodes in all files
 CALL imcoor (nimpf,lab0,impda,els_name,actio )

   !------  READ materials and sections
 DO i=1,nimpf                     !for each file
   fimpo = 70 + i                 !file to read from
   ! materials and sections
   IF( impda(5,i) ) THEN          !if mats and secs were stored
     READ (70+i)nm,ns   !number of materials and sections
     mn(:,:,i) = -1               !initializes relation
     sn(:,:,i) = -1
     CALL immtsc (fimpo,nm,ns,mn(:,:,i),sn(:,:,i))  !read mats and secs

   END IF
 END DO

 !------  READ element connectivities & internal variables
 DO i=1,nimpf                     !for each set
   fimpo = 70+i                   !file to read from
   lab1 = lab0(i)                 !label
   renum = impda(2,i)             !flag to change labels
   seque = impda(3,i)             !sequential or adding
   inter = impda(8,i)             !flag if internal variables are stored
   overw = impda(9,i)             !flag if and old set is overwrited
   snn = sn(:,:,i)     !section association ==> MAT_DBA
   !  re-read node labels of the node-set
   READ(fimpo) numpo              !number of nodes in the set
   IF( ASSOCIATED(nodset) )DEALLOCATE(nodset)
   ALLOCATE(nodset(numpo))
   READ(fimpo)(nodset(j),j=1,numpo)  !read node labels
   !  read element information
   READ (fimpo) itype                !read element type
   IF( itype > 30 )itype = (itype-MOD(itype,10))/10 !
   CALL elmdat ('IMPORT',j,els_name(i),itype)
 END DO
 CALL elsets ( )               !computes ESETS & ESET(:)

 !------  READ boundary conditions
 nd1 = ndofn                         !number of DOFS of the problems
 IF( rot_free )nd1 = nd1+1              !number of DOFs of the problem

 DO i=1,nimpf                    !for each set
   IF( impda(6,i) )THEN          !if boundary conditions are included
     fimpo = 70+i        !file to read from
     lab1 = lab0(i)      !label
     renum = impda(2,i)  !flag to change labels
     seque = impda(3,i)  !sequential or adding
     !  re-read node labels of the node-set
     READ(fimpo) numpo,idime,ieulr            !number of nodes in the set
     IF( ASSOCIATED(nodset) )DEALLOCATE(nodset)
     ALLOCATE(nodset(numpo))
     READ(fimpo)(nodset(j),j=1,numpo)
     ! read boundary conditions
     CALL impfix (idime,ieulr)
   END IF

   CLOSE(fimpo)          !close file
 END DO
 DEALLOCATE( lab0, impda, els_name, mn, sn, nodset )  !release memory of auxiliar arrays
 RETURN
 END SUBROUTINE import
