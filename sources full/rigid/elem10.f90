 SUBROUTINE elem10(TASK, nelms, elsnam, istop, flag1, flag2)

 ! Tasks Routine for element RIGID - HEAT

 USE ctrl_db, ONLY: ndime, ndofn, npoin, lumped
 USE outp_db, ONLY: sumat, iwrit
 USE ele10_db
 USE kinc_db
 USE npo_db
 USE surf_db


 IMPLICIT NONE
 CHARACTER(len=*),INTENT(IN):: TASK

 ! optional parameters

 LOGICAL, OPTIONAL :: flag1,flag2
 CHARACTER (len=*), OPTIONAL  :: elsnam
 INTEGER (kind=4), OPTIONAL  :: istop, nelms(:)

 CHARACTER(len=mnam) :: sname

 INTEGER (kind=4) nelem,nnode,nmast,ntype

 TYPE (ele10_set), POINTER :: elset, anter

 IF ( .NOT.ASSOCIATED (head) ) RETURN

 NULLIFY (anter)
 elset => head

 DO

   CALL comm10(1,nelem,ntype,nnode,nmast,sname,elset)


   SELECT CASE (TRIM(task))

   CASE ('GAUSSC')
     IF( elset%heat ) CALL gaus10(nelem, ntype, elset%ngaus, nnode, elset%lnods, coord,    &
                   elset%cartd, elset%dvolu, elset%shape, istop)

   CASE ('LUMASS')
     IF( ntype == 0 ) THEN   !defined by particles
       IF( elset%tmass /= 0d0 )THEN
       IF( nmast /= 0 ) CALL luma10p (ndime,nnode,elset%lnods,emass,sumat,elset%tmass)

       ELSE
         !nothing, Concentrated masses must be used
       END IF
     ELSE
       IF( nmast /= 0 ) CALL luma10 (nelem,ntype,ndime,nnode,elset%lnods,elset%matno,coord,  &
                                     emass,iwrit,sumat, elset%tmass, lumped,ymass,ifpre)
     END IF

   CASE ('NEW','NSTRA1','NSTRA2')
     CALL acvd10 (nnode,ndime,nelem,ndumn,elset%lnods,nmast,    &
                  ndofn,ifpre,elset%heat)

   CASE ('UPDLON')
     CALL updl10 (elset,oldlb)

   CASE ('RIGBDY')
     !compute center of mass and inertia tensor in principal directions
     IF( ntype == 0 ) THEN   !defined by particles
       CALL rigbdy(nmast,elset%lnods(:,1),nnode)
     ELSE
       IF( nmast /= 0 ) CALL rigb10(nelem,elset%lnods,nmast,nnode)
     END IF

   CASE ('DUMPIN')
     CALL dump10 (sname, nelem, nnode, elset%rbnod, nmast, ntype,   &
                  elset%lnods,elset%matno, elset%tmass,             &
                  elset%heat, elset%ngaus, elset%shape, elset%cartd, elset%dvolu)

   CASE ('WRTPOS')
     CALL mase10(nnode,nelem,elset%matno,elset%lnods,elset%sname,   &
                 elset%ntype, elset%nmast )

   CASE ('DELETE')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       CALL del_ele10 (head, tail, anter, elset)
       nelms(10) = nelms(10) - 1
       EXIT
     END IF
     IF ( ASSOCIATED (elset%next) ) anter => elset

   CASE ('SURFAC')
     IF( flag2 )EXIT
     flag2 = TRIM(elsnam) == TRIM(sname)
     IF (flag2) THEN
       IF( PRESENT(flag1)) THEN
         IF (ASSOCIATED(surfa%head)) CALL delete_seg(surfa%head,surfa%tail)
         NULLIFY(surfa%head,surfa%tail)          !Initialize empty list
         surfa%nelem = nmast
       ELSE
         ! get surface definition from the element set
         CALL surf10 (elset%lnods,elset%nelem,nnode,ndime,ntype)
       END IF
       EXIT
     END IF

   CASE ('TLMASS')      !compute (lumped) capacity vector
     IF( elset%heat )CALL tlma10 (nelem,ntype,ndime,nnode,elset%ngaus,elset%lnods,elset%dvolu,     &
                  elset%shape,elset%matno,tmass,iftmp)

   CASE ('TRESID')      !compute internal residual heats
     IF( elset%heat )CALL tres10(elset%ngaus,nnode,nelem,ndime,ntype,elset%matno,elset%lnods,   &
                 elset%shape,elset%cartd,elset%dvolu,tresi)

   END SELECT

   IF ( ASSOCIATED (elset%next) ) THEN
     elset => elset%next
   ELSE
     EXIT
   ENDIF

 END DO
 !     delete null elements sets after Nodal Updating
 IF( TRIM(task) == 'UPDLON')THEN
   NULLIFY (anter)
   elset => head
   DO
     IF(.NOT.ASSOCIATED(elset))EXIT    
     IF( elset%nelem == 0 )THEN
       CALL del_ele10 (head, tail, anter, elset)
       nelms(10) = nelms(10) - 1
     END IF
     IF (ASSOCIATED(elset))THEN
       anter => elset
       elset => elset%next
     ELSE IF(ASSOCIATED(anter))THEN
       elset => anter%next
     ELSE IF (ASSOCIATED(head)) THEN
       elset => head
     ELSE
       EXIT
     END IF
   END DO
 END IF

 RETURN
 END SUBROUTINE elem10
