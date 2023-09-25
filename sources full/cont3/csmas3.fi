 SUBROUTINE csmas3(lcseg,x,emass,nsegm,density)
 ! Compute nodal mass for contact surfaces if required
 IMPLICIT NONE
 !dummy arguments
 INTEGER (kind=4), INTENT(IN) :: nsegm,lcseg(:,:)
 REAL (kind=8), INTENT(IN) :: x(:,:),density
 REAL (kind=8), INTENT(IN OUT) :: emass(:,:)

 !local Variables
 INTEGER (kind=4) :: n,n1,n2,n3
 REAL (kind=8) :: t1(3),t2(3),t3(3),a2,mass

 INTERFACE
    INCLUDE 'vecpro.h'
    INCLUDE 'vecuni.h'
 END INTERFACE

 DO n=1,nsegm                   ! for each segment
   n1 = lcseg(1,n)              ! first node
   n2 = lcseg(2,n)              ! second node
   n3 = lcseg(3,n)              ! third node
   t1 = x(:,n2) - x(:,n1)       ! side 3
   t2 = x(:,n3) - x(:,n1)       ! side 2
   CALL vecpro(t1,t2,t3)        ! t3 = t1 x t2
   CALL vecuni(3,t3,a2)         ! a2 = twice the area
   mass = a2*density/6d0        ! nodal mass = Total_mass/3
   emass(1:3,n1) = emass(1:3,n1) + mass  !add to each node
   emass(1:3,n2) = emass(1:3,n2) + mass
   emass(1:3,n3) = emass(1:3,n3) + mass
 END DO

 RETURN
 END SUBROUTINE csmas3
