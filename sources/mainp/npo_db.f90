 MODULE npo_db
  ! nodal point information
  IMPLICIT NONE

  INTEGER (kind=4),POINTER :: &
            ifpre(:,:),   & !(ndofn,npoin) Associated equation to each DOF's
            iffix(:),     & !(npoin)       rotation constraint (shell elements)
            iftmp(:,:),   & !(ndoft,npoin) Associated equation to each Temp DOF
            label(:),     & !(npoin)       node labels
            loass(:),     & !(nload)       curve assigment for each load set
            heata(:),     & !(nheat)       curve assigment for each heat set
            oldlb(:),     & !(npoio)       node labels in previous strategy
            ifact(:)        !(npoin)       factors to compute coorb & coort

  REAL (kind=8),POINTER :: &
          acceg(:,:),     & !(nacce,ndime) ground accelerations
          coora(:,:),     & !(ndime,npoin) actual coordinates
          coorc(:,:),     & !(ndime,npoin) coordinates at previous step
          coord(:,:),     & !(ndime,npoin) original coordinates
          coors(:,:),     & !(ndime,npoin) coordinates at beginning of stage
          coorb(:,:),     & !(ndime,npoin) coordinates of bottom surface
          coort(:,:),     & !(ndime,npoin) coordinates of top surface
          emass(:,:),     & !(ndofn,npoin) nodal mass
          euler(:,:),     & !(neulr,npoin) local coordinate system (actual)
          eule0(:,:),     & !(1-3,npoin) local coordinate system (initial)
          fcont(:,:),     & !(ndime,npoin) contact forces
          force(:,:),     & !(neq,nload+1) external forces
          loadv(:,:,:),   & !(ndofn,npoin,nload) external forces
          resid(:,:),     & !(ndofn,npoin) internal forces
          psic(:,:),      & !(2,npoin) psi functions (beginning of stage)
          psia(:,:),      & !(2,npoin) psi functions (present)
          heati(:,:),     & !(neqt,nheat+1) external heating
          tresi(:),       & !(neqt) heat transfer residual
          dtemp(:,:),     & !(ndoft,npoin) temperature increments
          tempe(:,:),     & !(ndoft,npoin) temperature
          velnp(:,:)        !(ndofn,npoin) nodal velocities (assigned at outdyn)

  LOGICAL ,POINTER :: naeul(:) !(npoin) if local system is time dependent

   REAL (kind=8),POINTER ::&
        acelr(:), & !(neq) accelerations
        ddisp(:), & !(neq) incremental displacements
        veloc(:), & !(neq) velocities
        trate(:), & !(neqt) temperature rate
        tmass(:), & !(neqt) equivalent DOF Thermal Capacity for nodes
        ymass(:)    !(neq) equivalent DOF mass

 CONTAINS

   SUBROUTINE dump_npo

     !   dumps database of control parameters

     USE ctrl_db, ONLY: ndime, ndofn, neulr, nload, npoin, nacce, neq, numct, ndimc, &
                        bottom, top, itemp, ndoft, nheat, therm, neqt
     USE esets_db,ONLY: nelms, rot_free

     IMPLICIT NONE

     INTEGER (kind=4) :: i,j,k
     LOGICAL :: psi

     WRITE (50,ERR=9999) ((ifpre(i,k),i=1,ndofn),k=1,npoin)

     IF( rot_free )THEN
       WRITE(50,ERR=9999)( iffix(k),k=1,npoin)
     ELSE
       WRITE(50,ERR=9999) iffix(1)
     END IF

     WRITE(50,ERR=9999) (label(i),i=1,npoin)
     WRITE(50,ERR=9999) ((coora(i,k),i=1,ndime),k=1,npoin), &
                        ((coorc(i,k),i=1,ndime),k=1,npoin), &
                        ((coord(i,k),i=1,ndime),k=1,npoin), &
                        ((coors(i,k),i=1,ndime),k=1,npoin), &
                        ((emass(i,k),i=1,ndofn),k=1,npoin), &
                        ((resid(i,k),i=1,ndofn),k=1,npoin)
     WRITE(50,ERR=9999) ((acceg(i,j),i=1,nacce),j=1,ndime)
     IF(nload > 0)WRITE(50,ERR=9999) ((force(i,k),i=1,neq+1),k=1,nload+1)

      IF( itemp ) THEN
        WRITE(50,ERR=9999) ((tempe(k,i),k=1,ndoft),i=1,npoin)
        WRITE(50,ERR=9999) ((dtemp(k,i),k=1,ndoft),i=1,npoin)
        WRITE(50,ERR=9999) ((iftmp(k,i),k=1,ndoft),i=1,npoin)
        IF( therm )THEN
          IF( nheat > 0)  WRITE(50,ERR=9999) (heata(k),k=1,nheat)
          WRITE(50,ERR=9999) ((heati(k,i),k=1,neqt),i=1,nheat+1)
          WRITE(50,ERR=9999) (tresi(k),k=1,neqt)
        END IF
      END IF

      IF(neulr > 0)THEN
        WRITE(50,ERR=9999)(eule0(:,k),k=1,npoin)
        WRITE(50,ERR=9999)((euler(i,k),i=1,neulr),k=1,npoin),(naeul(k),k=1,npoin)
        psi = ASSOCIATED(psia)
        WRITE(50,ERR=9999)psi
        IF( psi )THEN
          WRITE(50,ERR=9999)((psic(i,j),i=1,2),j=1,npoin),((psia(i,j),i=1,2),j=1,npoin)
        END IF
      ELSE
        WRITE(50,ERR=9999) euler(1,1)
      END IF

      IF(numct > 0)THEN
        WRITE(50,ERR=9999) ((fcont(i,k),i=1,ndimc),k=1,npoin)
      ELSE
        WRITE(50,ERR=9999) fcont(1,1)
      END IF

      IF( nload > 0)THEN
        WRITE(50,ERR=9999)loass
        WRITE(50,ERR=9999)(((loadv(i,j,k),i=1,ndofn),j=1,npoin),k=1,nload)
      END IF

      WRITE(50,ERR=9999) ((velnp(i,k),i=1,ndofn),k=1,npoin)

      WRITE(50,ERR=9999) (acelr(i),i=1,neq),(ddisp(i),i=1,neq),                  &
                         (veloc(i),i=1,neq),(ymass(i),i=1,neq)
      IF( therm ) WRITE(50,ERR=9999) (trate(i),i=1,neqt),(tmass(i),i=1,neqt)

   RETURN
   9999 CALL runen2('')
   END SUBROUTINE dump_npo

   SUBROUTINE rest_npo

     !   restores database

     USE ctrl_db, ONLY: ndime, ndofn, neulr, nload, npoin, neq, nacce, numct, &
                        ndimc, bottom, top, therm, itemp, ndoft, nheat, neqt
     USE esets_db

     INTEGER (kind=4) :: i,j,k
     LOGICAL :: psi

     ALLOCATE( ifpre(ndofn,npoin), label(npoin),   &
               velnp(ndofn,npoin), oldlb(npoin))
     READ (51) ((ifpre(i,k),i=1,ndofn),k=1,npoin)

     IF( rot_free )THEN
       ALLOCATE( iffix(npoin) )
       READ (51)( iffix(k),k=1,npoin)
     ELSE
       ALLOCATE( iffix(1) )
       READ (51) iffix(1)
     END IF

     READ (51) (label(i),i=1,npoin)
     DO i=1,npoin
       oldlb(i) = label(i)
     END DO
     ALLOCATE( coora(ndime,npoin), coorc(ndime,npoin),           &
             coord(ndime,npoin), coors(ndime,npoin),             &
             emass(ndofn,npoin),                                 &
             resid(ndofn,npoin))
     READ (51) ((coora(i,k),i=1,ndime),k=1,npoin), &
               ((coorc(i,k),i=1,ndime),k=1,npoin), &
               ((coord(i,k),i=1,ndime),k=1,npoin), &
               ((coors(i,k),i=1,ndime),k=1,npoin), &
               ((emass(i,k),i=1,ndofn),k=1,npoin), &
               ((resid(i,k),i=1,ndofn),k=1,npoin)

       ALLOCATE( acceg(nacce,ndime))
       READ(51) ((acceg(i,k),i=1,nacce),k=1,ndime)
       IF(nload > 0)THEN
         ALLOCATE(force(neq+1,nload+1))
         READ(51) ((force(i,k),i=1,neq+1),k=1,nload+1)
       END IF

      IF( itemp ) THEN
        ALLOCATE( tempe(ndoft,npoin), dtemp(ndoft,npoin), iftmp(ndoft,npoin) )
        READ(51) ((tempe(k,i),k=1,ndoft),i=1,npoin)
        READ(51) ((dtemp(k,i),k=1,2),i=1,npoin)
        READ(51) ((iftmp(k,i),k=1,ndoft),i=1,npoin)
        IF( therm )THEN
          ALLOCATE( heata(nheat))
          IF( nheat > 0) READ(51) (heata(k),k=1,nheat)
          ALLOCATE( heati(neqt,nheat+1), tresi(neqt) )
          READ(51) ((heati(k,i),k=1,neqt),i=1,nheat+1)
          READ(51) (tresi(k),k=1,neqt)
        END IF
      ELSE
        ALLOCATE( tempe(1,1), dtemp(1,1) )
      END IF

     IF(neulr > 0)THEN
       ALLOCATE(eule0(2*ndime-3,npoin))
       ALLOCATE(euler(neulr,npoin),naeul(npoin))
       READ (51)(eule0(:,k),k=1,npoin)
       READ (51)((euler(i,k),i=1,neulr),k=1,npoin),(naeul(k),k=1,npoin)
       READ(51)psi
       IF( psi )THEN
         ALLOCATE(psic(2,npoin),psia(2,npoin))
         READ(51)((psic(i,j),i=1,2),j=1,npoin),((psia(i,j),i=1,2),j=1,npoin)
       END IF
     ELSE
       ALLOCATE( euler (1,1))
       READ (51) euler(1,1)
     END IF

     IF(numct > 0)THEN
       ALLOCATE ( fcont(ndimc,npoin) )
       READ (51) ((fcont(i,k),i=1,ndimc),k=1,npoin)
     ELSE
       ALLOCATE ( fcont(1,1) )
       READ (51) fcont(1,1)
     END IF

     IF( nload > 0)THEN
       ALLOCATE(loass(nload))
       READ(51)loass
       ALLOCATE( loadv(ndofn,npoin,nload) )
       READ (51)(((loadv(i,j,k),i=1,ndofn),j=1,npoin),k=1,nload)
     ELSE
       NULLIFY(force,loass,loadv)
     END IF

     READ (51) ((velnp(i,k),i=1,ndofn),k=1,npoin)

     ALLOCATE( acelr(neq),ddisp(neq),veloc(neq),ymass(neq))
     READ (51) (acelr(i),i=1,neq),(ddisp(i),i=1,neq),   &
               (veloc(i),i=1,neq),(ymass(i),i=1,neq)
     IF( therm )THEN
       ALLOCATE( trate(neqt), tmass(neqt))
       READ(51) (trate(i),i=1,neqt),(tmass(i),i=1,neqt)
     END IF

   END SUBROUTINE rest_npo
 END MODULE npo_db
