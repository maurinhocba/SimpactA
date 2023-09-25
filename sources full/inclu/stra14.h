 SUBROUTINE stra14(a,b,x,t3,t,h,sides,is,stran,lb)

 !     Compute first and second fundamental forms for element CST-BST (TLF)
 !USE ele14_db, ONLY : kk,hh
 IMPLICIT NONE

 REAL (kind=8), INTENT(IN) :: x(3,6),    &  !nodal cooordinates
                              b(3,0:3),  &  !cartesian derivatives local x1
                              a(3,0:3)      !cartesian derivatives local x2
 REAL (kind=8), INTENT(OUT) :: t(3,2,0:3), & !derivatives of the element configuration
                               t3(3),      & !element normal
                               h(3,3)        !h vectors
 LOGICAL, INTENT(IN) :: sides(3),  & !True = side elements exists
                        is(3)        !True = clamped side
 REAL (kind=8), INTENT(OUT), OPTIONAL :: stran(6), & !1st & 2nd fund. forms
                                         lb            !thickness ratio
 END SUBROUTINE stra14
