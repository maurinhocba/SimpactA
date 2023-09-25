      SUBROUTINE aplfol (ndime,coora,resid,lcurv,dtime,ttime, &
     &                   ntype,fltype,flparm,headf,factor)

      !     applies follower load

      USE curv_db
      USE loa_db
      IMPLICIT NONE
      CHARACTER (len=*) :: fltype
      INTEGER (kind=4) :: lcurv,ndime,ntype
      REAL (kind=8) :: dtime,ttime,resid(:,:),coora(:,:),factor
      TYPE (flpar) :: flparm
      TYPE (foll_seg),POINTER :: headf

      END SUBROUTINE aplfol
