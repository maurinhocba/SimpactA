 SUBROUTINE stgp06(ngaus,s,e,d)
 !*****************************************************************************
 !
 !*****evaluates total resultant stresses for shell element
 !     for linear elastic isotropic material
 !
 !****************************************************************************
 IMPLICIT NONE

 !                        routine parameters

 INTEGER (kind=4), INTENT(IN) :: ngaus

 REAL (kind=8), INTENT(IN) :: e(:,:),d(:)
 REAL (kind=8), INTENT(OUT) :: s(:,:)

 !                        local variables
 INTEGER (kind=4) :: g  !Gauss point Index

 DO g = 1,ngaus
   !                     compute effective stresses
 !             1  2  3 16 17 18        35 36 37 38
 !             2  4  5 19 20 21        39 40 41 42
 !             3  5  6 22 23 24        43 44 45 46
 !            16 19 22  7  8  9        47 48 49 50
 !            17 20 23  8 10 11        51 52 53 54
 !            18 21 24  9 11 12        55 56 57 58
 !                              13 14              62 63
 !                              14 15              64 65
 !            35 39 43 47 51 55        25 26 27 28
 !            36 40 44 48 52 56        26 29 30 31
 !            37 41 45 49 53 57        27 30 32 33
 !            38 42 46 50 54 58        28 31 33 34
 !                              62 64              59 60
 !                              63 65              60 61
   !  membrane forces
   s(1,g) = DOT_PRODUCT(d((/ 1, 2, 3,16,17,18,35,36,37,38/)),e((/1:6,9:12/),g))
   s(2,g) = DOT_PRODUCT(d((/ 2, 4, 5,19,20,21,39,40,41,42/)),e((/1:6,9:12/),g))
   s(3,g) = DOT_PRODUCT(d((/ 3, 5, 6,22,23,24,43,44,45,46/)),e((/1:6,9:12/),g))
   ! bending moments
   s(4,g) = DOT_PRODUCT(d((/16,19,22, 7, 8, 9,47,48,49,50/)),e((/1:6,9:12/),g))
   s(5,g) = DOT_PRODUCT(d((/17,20,23, 8,10,11,51,52,53,54/)),e((/1:6,9:12/),g))
   s(6,g) = DOT_PRODUCT(d((/18,21,24, 9,11,12,55,56,57,58/)),e((/1:6,9:12/),g))
   ! additional moments
   s( 9,g)= DOT_PRODUCT(d((/35,39,43,47,51,55,25,26,27,28/)),e((/1:6,9:12/),g))
   s(10,g)= DOT_PRODUCT(d((/36,40,44,48,52,56,26,29,30,31/)),e((/1:6,9:12/),g))
   s(11,g)= DOT_PRODUCT(d((/37,41,45,49,53,57,27,30,32,33/)),e((/1:6,9:12/),g))
   s(12,g)= DOT_PRODUCT(d((/38,42,46,50,54,58,28,31,33,34/)),e((/1:6,9:12/),g))
   ! shear forces
   s(7,g) = DOT_PRODUCT(d((/13,14,62,63/)),e((/7,8,13,14/),g))
   s(8,g) = DOT_PRODUCT(d((/14,15,64,65/)),e((/7,8,13,14/),g))
   ! additional shear forces
   s(13,g) = DOT_PRODUCT(d((/62,64,59,60/)),e((/7,8,13,14/),g))
   s(14,g) = DOT_PRODUCT(d((/63,65,60,61/)),e((/7,8,13,14/),g))
 END DO
 RETURN
 END SUBROUTINE stgp06
