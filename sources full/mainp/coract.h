 SUBROUTINE coract(lambd,fixvl,ddisp,coora,euler)

 !     UPDATES CONFIGURATION AND NODAL SYSTEMS

 !USE ctrl_db, ONLY : ndime,npoin,ndofn,neulr
 !USE kinc_db, ONLY : nn,nn1,nm,nescv,ndepd,naris,&
 !                    nndpd,distd,npsdf,nesdf,ftsdf
 !USE npo_db, ONLY : naeul,ifpre
 !USE rfdp_db, ONLY : nrfdp
 IMPLICIT NONE
 REAL (kind=8),INTENT(IN) :: ddisp(:),lambd,fixvl(:)
 REAL (kind=8),INTENT(IN OUT) :: coora(:,:)
 REAL (kind=8), POINTER :: euler(:,:)
 END SUBROUTINE coract
