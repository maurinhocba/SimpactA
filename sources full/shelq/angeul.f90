 SUBROUTINE angeul( romtx, angls, rad )
 !*************************************************************************
 !
 !     transformation from rotation matrix to euler's angles
 !
 !              xg = r xn  : ang glob -> nodl
 !
 !     angls(1) : rotation about z axe to take x axis in X-Y plane
 !     angls(2) : rotation about x axe to make z and Z coincide
 !     angls(3) : rotation about z axe to make x-y coincide with X-Y plane
 !
 !*************************************************************************
 IMPLICIT NONE
 !dummy arguments
 REAL (kind=8), INTENT(IN) :: romtx(3,3)  !rotation matrix (t1,t2,t3)
 REAL (kind=8), INTENT(OUT) :: angls(3)   !Euler angles
 LOGICAL, INTENT(IN), OPTIONAL :: rad     !if present angles returned in rads
 !local variables
 REAL (kind=8)  cosda,cosdg,sinda,sindg
 REAL (kind=8), PARAMETER :: pid2  = 1.5707963267948966192313216916398, &   ! pid2 = ATAN(1d0)*2d0
                             pi    = 3.1415926535897932384626433832795, &   ! pi   = 2d0*pid2
                             facto = 57.295779513082320876798154814105, &   ! facto= 180d0/pi
                             cuzer =1.0d-12


 IF ( romtx(3,3) >= 1d0 ) THEN          !full positive proyection over Z axis
   angls(2) = 0d0                           !second angle is null
 ELSE IF ( romtx(3,3) <= -1d0) THEN     !full negative proyection over Z axis
   angls(2) = pi                            !second angle is Pi
 ELSE                                   ! an intermediate angle exist
   angls(2) = ACOS( romtx(3,3) )            !compute angle
 END IF
 IF (angls(2) == 0d0) THEN              !IF t3 is along Z axis
   angls(1) = 0d0                           !first angle is null
 ELSE
   IF ( ABS(romtx(2,3)) < cuzer ) THEN
     IF (romtx(1,3) > 0d0) THEN
        angls(1) = pid2
     ELSE
        angls(1) =-pid2
     END IF
   ELSE
     angls(1) = ATAN( -romtx(1,3)/romtx(2,3) )
     IF (romtx(2,3) > 0d0 ) angls(1) = angls(1) + pi  !  cosine of beta negative
   END IF
 END IF
 cosda = COS(angls(1))
 sinda = SIN(angls(1))
 cosdg = romtx(1,1)*cosda + romtx(2,1)*sinda
 sindg =-romtx(1,2)*cosda - romtx(2,2)*sinda
 IF (ABS(cosdg) < cuzer) THEN
   IF (sindg > 0d0) THEN
     angls(3) = pid2
   ELSE
     angls(3) =-pid2
   END IF
 ELSE
   angls(3) = ATAN( sindg/cosdg )
   IF (cosdg < 0.0) angls(3) = angls(3) + pi      ! cosine negative
 END IF

 ! if code rad is not used in calling
 IF (.NOT.PRESENT (rad))  angls = angls*facto     ! change radians to degrees

 RETURN
 END SUBROUTINE angeul
