SUBROUTINE slmass ()

  ! *** calculates lumped mass for each element & assembles global mass matrix
  ! to be used in STATIC approach

  USE cms_db
  USE ctrl_db, ONLY : ndofn, npoin, neq
  USE npo_db,  ONLY : ymass,ifpre,emass
  USE kinc_db, ONLY : npsdf,nesdf,ftsdf
  USE sms_db,  ONLY : selective_mass_scaling, fibre_mass
  USE static_db, ONLY : smass

  IMPLICIT NONE

  INTEGER (kind=4) :: i,ipoin,chnode
  REAL    (kind=8) :: xcmas(ndofn),factor

  INTERFACE
    INCLUDE 'elemnt.h'
  END INTERFACE

  smass = 0d0             !initializes modified mass matrix
  CALL elemnt ('SLUMAS')  !compute element modified masses

  ! compute average factor for mass scaling
  i = 0
  factor = 0d0
  DO ipoin=1,npoin
    IF( smass(1,ipoin) > 0 )THEN
      i = i+1
      factor = factor + smass(1,ipoin) / emass(1,ipoin)
    END IF
  END DO
  factor = factor/i  !average factor

  !add concentrated masses
  IF (nconm > 0) THEN
    DO i = 1,nconm
      ipoin = nodcms(i)
      ipoin = chnode(ipoin)
      xcmas(1:ndofn) = cmass(1:ndofn,i)
      smass(1:ndofn,ipoin) = smass(1:ndofn,ipoin) + xcmas*factor
    END DO
  END IF

  !add other masses (contact)
  DO ipoin=1,npoin
    IF( smass(1,ipoin) == 0d0 )THEN
      IF( emass(1,ipoin) /= 0d0 ) smass(:,ipoin) = emass(:,ipoin)*factor
    END IF
  END DO
  !Assembles mass for active DOFs
  ymass = 0d0
  CALL ensmal(ndofn*npoin,neq,ifpre(1,1),smass(1,1),ymass(1),npsdf(1),nesdf(1),ftsdf(1))
  IF( selective_mass_scaling )CALL fibre_mass(smass)
  DO i=1,neq
    ymass(i) = 1d0/ymass(i)
  END DO
RETURN
9999 CALL runen2('')
END SUBROUTINE slmass
