 SUBROUTINE elem20(TASK, elsnam, dtime, ttime, istop, ivect, &
                   lnod, flag1, flag2)

 !master routine for element 20 (TLF) 2-D solid triangular element

 USE ctrl_db, ONLY: ndofn, npoin, ntype
 USE outp_db, ONLY: sumat, iwrit
 USE ele20_db
 USE npo_db
 USE sms_db, ONLY : selective_mass_scaling, sms_ns , sms_name, sms_thl, sms_alp
 IMPLICIT NONE

   !--- Dummy variables
   CHARACTER(len=*),INTENT(IN):: TASK
     ! optional parameters
   LOGICAL, OPTIONAL :: flag1,flag2
   CHARACTER (len=*), OPTIONAL :: elsnam
   INTEGER (kind=4), OPTIONAL :: istop, ivect(:)
   INTEGER (kind=4), POINTER, OPTIONAL :: lnod(:,:)
   REAL (kind=8), OPTIONAL :: dtime,ttime
   !--- Local variables
   CHARACTER (len=mnam) :: sname
   INTEGER (kind=4) :: nreqs, narch, nelem, i
   TYPE (ele20_set), POINTER :: elset, anter
   LOGICAL :: sms,esta
   REAL (kind=8) :: thl,alp

   IF (.NOT.ASSOCIATED(head)) RETURN  !no element sets

   NULLIFY(anter)      !nullify pointer to previous set
   elset => head        !point to first set of elements

   DO                   !loop over the element sets
     ! recover control variables for present element set
     CALL comm20(1,nelem,nreqs,narch,sname,elset)

     SELECT CASE(TRIM(task))   !according to the requested task

     CASE('GAUSSC')      !compute initial constants
       CALL gaus20(elset%head,coord,istop,elset%gauss,elset%angdf,ntype)

     !CASE('LOADPL')      !compute gravity load vector
     !  CALL load20(igrav,loadv(:,:,iload),gvect,gravy,elset%head,ntype,coord)

     CASE('CLOSEF')      !close output file
       CALL close1(nreqs,narch)

     CASE('LUMASS')      !compute (lumped) mass vector
       CALL luma20(elset%head, emass, sumat, ntype, coord)

     CASE('NEW','NSTRA1','NSTRA2')  !release DOFs and other tasks
       CALL acvd20(ifpre, elset%head, elset%lside, elset%swapc, elset%eside)
       CALL dist20(elset%head, elset%nelem, elset%eside, coord, .FALSE.)

     CASE('UPDLON')      !update local node numbers
       CALL updl20(elset%head, oldlb, elset%lside, elset%origl)

     CASE('OUTDYN')      !output variables for post-processing
       IF (flag1 .OR. flag2) CALL outd20(flag1,flag2,iwrit,elset%head,nreqs,narch,   &
                                         elset%ngrqs,ttime)

     CASE('RESVPL')      !compute internal nodal forces
       CALL resv20(nelem,elset%head,ntype,coora,resid,istop,ttime,elset%angdf,coorc,elset%eulrf)

     CASE('DUMPIN')      !dumps variables for restart
       CALL dump20(elset)

     CASE('WRTPOS')      !writes geometric variables for post-process
       CALL mase20(ntype,nreqs,nelem,elset%head,elset%ngrqs,narch,elset%angdf,elset%sname)
       elset%narch = narch              !keep output file

     CASE('INCDLT')      !compute critical time increment
       sms = .FALSE.
       thl = 0d0
       IF( selective_mass_scaling )THEN
         DO i=1,sms_ns
           IF(TRIM(sms_name(i)) == TRIM(sname) )THEN
             alp = sms_alp(i)
             thl = sms_thl(i)
             sms = .TRUE.
             EXIT
           END IF
         END DO
       END IF
       CALL delt20( nelem, elset%head, dtime, coora, sms, alp)

     CASE('SEARCH')      !search if set named SNAME exists
       IF (flag2) EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) EXIT

     CASE('SURFAC')
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! get surface definition from the element set
         CALL surf20(elset%head,elset%nelem)
         EXIT
       END IF

     CASE('BOUNDA')
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! get boundary definition from the element set for remesh
         CALL boun20(elset%head,elset%nelem)
         EXIT
       END IF

     CASE('SLNODS')   !extract connectivities into array
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         esta = PRESENT(flag1)
         IF( esta ) esta = flag1
         CALL slnods20(elset%head,lnod,nelem,esta)
         EXIT
       END IF

     CASE('SMOOTH')      !performs smoothing for remeshing and transfer
       IF (flag2) EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! smooth gaussian variables and writes into STRND
         CALL smth20(npoin,strnd,elset%head,elset%eulrf)
         EXIT
       END IF

     CASE('NODSET')      !compute nodes set from element set
       IF (flag2) EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! extracts nodes used in the discretization to nodset
         CALL nods20(npoin, numpo, elset%head, nodset, label)
         EXIT
       END IF

     CASE ('SECDAT')
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         !  extract SECTIONs used in the set
         CALL secd20 (elset%head,elset%nelem,ivect)
         EXIT
       END IF

     CASE ('EXPORT')
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! export set data
         CALL expo20 (elset,flag1,istop)
         EXIT
       END IF

     CASE('INTVAR')
       IF (flag2) EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ALLOCATE(sg_old(7,1,elset%nelem)) !7:number of transfer variables by element
         ! get element information to transfer
         CALL getelv20(elset%head,sg_old,elset%eulrf)
         EXIT
       END IF

     CASE('REALLO') !only for the remeshed set
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         ! reallocates the element set after remeshing
         ! exchange mesh - deallocate old mesh and put new mesh,
         ! interpolate all the Gaussian variables from nodes to Gauss points
         elset%gauss = .FALSE.
         elset%lside = .FALSE.
         elset%plstr = -1
         CALL real20(elset%nelem,elset%head,elset%tail,iwrit,elset%eulrf)
         elset%origl = .TRUE.   !connectivities changed to original labels
         ! restore the neighbour element patch array
         DEALLOCATE(elset%eside)
         ALLOCATE(elset%eside(3,elset%nelem))
         elset%eside = 0     !reinitializes
       END IF


     !CASE('INIGAU')      !read initial stresses or internal variables
     !  IF( flag2 )EXIT
     !  flag2 = TRIM(elsnam) == TRIM(sname)
     !  IF (flag2) THEN
     !    ! read initial internal variables for the element set
     !    CALL inig20 (elset%head,nelem,iwrit,elset%plstr)
     !    EXIT
     !  END IF

     CASE('ELMDIS') !computes element distorsion for remeshing purposes
       IF( flag2 )EXIT
       flag2 = TRIM(elsnam) == TRIM(sname)
       IF (flag2) THEN
         CALL dist20(elset%head, elset%nelem, elset%eside, coora, .TRUE.)
         EXIT
       END IF

     END SELECT

     IF ( ASSOCIATED(elset%next)) THEN   !more sets to process
       elset => elset%next                 !point to next set
     ELSE
       EXIT
     END IF

   END DO

 RETURN

 END SUBROUTINE elem20
