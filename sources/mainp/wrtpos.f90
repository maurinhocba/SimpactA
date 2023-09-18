 SUBROUTINE wrtpos(flag)
 !
 ! write initial data for post-process
 ! updates stage coordinates for a new strategy
 !
 USE ctrl_db, ONLY: ndime,npoin,ndofn,nrotd,numct,text,ndoft,itemp
 USE c_input, ONLY: openfi
 USE npo_db, ONLY: label,coord,coors
 USE outp_db
 USE mat_dba, ONLY : nusect,wrt_mat_db
 USE flc_db, ONLY : nflc
 USE gvar_db, ONLY : static
 IMPLICIT NONE

   ! Dummy variables
   INTEGER(kind=4):: flag
   ! Local variables
   INTEGER (kind=4) :: i , idyna

   INTERFACE
     INCLUDE 'contac.h'
     INCLUDE 'elemnt.h'
   END INTERFACE

   IF (flag /= 1) THEN  !Open files and write initial data for post-processing
     CALL openfi(16)    ! open post-process files
     CALL openfi(17)
     ! open file for especific output (for HISTORY or CURVES)
     IF (nreqd > 0) THEN              !displacement
       CALL openfi(11)
       WRITE(11,ERR=9999)postype, nreqd,ndofn,(label(nprqd(i)),i=1,nreqd)
     END IF
     IF (nreql > 0) THEN              !equivalent nodal forces
       CALL openfi(12)
       WRITE(12,ERR=9999)postype, nreql,ndofn,(label(nprql(i)),i=1,nreql)
     END IF
     IF (iener == 1 .AND. .NOT.static ) THEN             !energies
       CALL openfi(15)
       WRITE(15,ERR=9999)postype, 11   !number of values to print per step
     END IF
     IF (nreqv > 0 .AND. .NOT.static ) THEN              !velocities
       CALL openfi(19)
       WRITE(19,ERR=9999)postype, nreqv,nrotd,(label(nprqv(i)),i=1,nreqv)
     END IF
     IF (nreqa > 0 .AND. .NOT.static ) THEN              !accelerations
       CALL openfi(20)
       WRITE(20,ERR=9999)postype, nreqa,nrotd,(label(nprqa(i)),i=1,nreqa)
     END IF
     IF (nreqt > 0) THEN              !temperatures
       CALL openfi(70)
       WRITE(70,ERR=9999)postype, nreqt,ndoft,(label(nprqt(i)),i=1,nreqt)
     END IF
      !writes information for post-processing associated to follower loads
     CALL wrtfl0()
     CALL actcd()  !updates stage coordinates
     ! write some control parameter
     i = 0
     IF( itemp ) i = ndoft
     idyna = 1                 !default analysis type
     IF( static ) idyna = 0    !static analysis type
     WRITE(17,ERR=9999) ndime,ndofn,i !i = number of termal DOFs
     WRITE(17,ERR=9999) npoin,idyna,text    !1=dynamic analysis, 0=Static analysis TEXT = problem title
     DO i=1,npoin
       WRITE(17,ERR=9999) label(i),coord(1:ndime,i),coors(1:ndime,i)
     END DO

     WRITE(17,ERR=9999) nusect
     CALL wrt_mat_db( )
     CALL elemnt ('WRTPOS') !write element data (mesh connectivities)
     WRITE(17,ERR=9999) 0,0,'                              '       !to indicate end of elements sets

     IF (numct > 0) THEN                 !contact
       IF (nreqc > 0) THEN               !contact forces at selected points
         CALL openfi(14)
         WRITE(14,ERR=9999) nreqc,ndime,(label(nprqc(i)),i=1,nreqc)
       END IF
       ! writing rigid (contact) surfaces for postprocess
       CALL contac ('WRTSUR',flag)
     END IF

   ELSE !IF (flag == 1) THEN   !Open files to add data for post-processing
 ! ***
     CALL openfi(16,flag=flag)      ! open post-process files
     CALL openfi(17,flag=flag)

     IF (nreqd > 0) CALL openfi(11,flag=flag)
     IF (nreql > 0) CALL openfi(12,flag=flag)
     IF (iener == 1 .AND. .NOT.static ) CALL openfi(15,flag=flag)
     IF (nreqv > 0 .AND. .NOT.static ) CALL openfi(19,flag=flag)
     IF (nreqa > 0 .AND. .NOT.static ) CALL openfi(20,flag=flag)
     IF (nreqt > 0) CALL openfi(70,flag=flag)
     IF (numct > 0) THEN
       IF (nreqc > 0) CALL openfi(14,flag=flag)
       CALL contac ('WRTSUR',iwrit=flag)
     END IF
 !        IF( thickc )CALL openfi(46,flag=flag)
 !        IF(nreqp > 0) CALL openfi(150,flag=flag)
 !        IF(nreqt > 0) CALL openfi(151,flag=flag)

   END IF

 RETURN
 9999 CALL runen2('')
 END SUBROUTINE wrtpos
