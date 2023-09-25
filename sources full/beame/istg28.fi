 SUBROUTINE istg28(ng,nstre,stres,ds,sedef,c)
 !**********************************************************************
 !
 !      incremental stresses at gauss points for beams (model 2)
 !
 !**********************************************************************
 IMPLICIT NONE
 INTEGER (kind=4) ng,nstre
 REAL    (kind=8) stres(nstre,ng),ds(nstre,ng), &
                  sedef(ng),c(:)
 INTEGER (kind=4) g

 DO g=1,ng
   sedef(g) = sedef(g) + ds(1,g)
   stres(1,g) =  +c( 5)*ds(1,g) +c(11)*ds(2,g) +c(12)*ds(3,g) &
                 +c(13)*ds(4,g) +c(14)*ds(5,g) +c(15)*ds(6,g)
   stres(2,g) =  +c(11)*ds(1,g) +c( 6)*ds(2,g) +c(16)*ds(3,g) &
                 +c(17)*ds(4,g) +c(18)*ds(5,g) +c(19)*ds(6,g)
   stres(3,g) =  +c(12)*ds(1,g) +c(16)*ds(2,g) +c( 7)*ds(3,g) &
                 +c(20)*ds(4,g) +c(21)*ds(5,g) +c(22)*ds(6,g)
   stres(4,g) =  +c(13)*ds(1,g) +c(17)*ds(2,g) +c(20)*ds(3,g) &
                 +c( 8)*ds(4,g) +c(23)*ds(5,g) +c(24)*ds(6,g)
   stres(5,g) =  +c(14)*ds(1,g) +c(18)*ds(2,g) +c(21)*ds(3,g) &
                 +c(23)*ds(4,g) +c( 9)*ds(5,g) +c(25)*ds(6,g)
   stres(6,g) =  +c(15)*ds(1,g) +c(19)*ds(2,g) +c(22)*ds(3,g) &
                 +c(24)*ds(4,g) +c(25)*ds(5,g) +c(10)*ds(6,g)
 END DO

 END SUBROUTINE istg28
