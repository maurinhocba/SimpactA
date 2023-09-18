  SUBROUTINE cohill(ppr,yieldf,assoc)
  ! computes coefficients for Hill 48 yield criteria
  ! from strength ratios and Lankford parameters
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN OUT) :: ppr(:)
  INTEGER(kind=4), INTENT(IN) :: yieldf,assoc
  ! PPR: (IN)    1-7: not used
  !             8-11:S90,Biuni,SXY,S45
  !            12-15:rM,r0,r45,r90
  !      (OUT) 16-21: Hill coefficients for yield function
  !            22-27: Hill coefficients for potential function

  REAL (kind=8) :: rm,r0,r45,r90,s90,s45,biuni,sxy,f,g,h,l,m,n
  REAL (kind=8) :: lb,ss90,ss45,sb,aux,aux1,mh

  LOGICAL ::  s,r

  SELECT CASE (yieldf)
  CASE (1) !elastic
    RETURN
  CASE (2) ! Mises coefficients
    ppr(16:18) = 0.5d0
    ppr(19:21) = 1.5d0
    ppr(22:27) = ppr(16:21)    !associative plasticity
  CASE (3) !Hill 48 anisotropic (PLANAR) plasticity
    rm = ppr(12)                     !average lankford ratio
    r = ANY(ppr(12:15) > 0d0 )       !Lankford ratios information exists
    s = ppr(8) /= 0d0                !Strength ratios information exists

    IF( s )THEN  ! relative strength data exists
      s90 = ppr(8)
      s45 = ppr(11)
      sxy = ppr(10)
      biuni = ppr(9)
      ! compute missing parameters
      IF( r .AND. s45 /= 0d0 .AND. biuni == 0d0) THEN
        ! IF s45 and no biaxial compute Lambda parameter
        IF( rm == 0d0 )rm = (ppr(13)+ppr(14)*2+ppr(15))/4d0
        CALL coef_l(rm,s90,s45,lb)
        rm = ppr(12)                     !average lankford ratio
        IF(lb > 0.7 ) WRITE(*,"('Hill coeficients may be poorly computed, CHECK THEM')") 
        sxy = SQRT(3d0)/SQRT( 4d0/s45**2+lb-1d0-1d0/s90**2)
      ELSE           ! biuni or SXY
        IF( biuni /= 0d0 )THEN
          lb = s90*(1d0+1d0/s90**2-1d0/biuni**2)
        ELSE
          lb = 1d0               !no biuni make Lambda = 1
        END IF
        IF(sxy == 0d0 )THEN              !if NO exists SXY
          IF(s45 /= 0d0)THEN
            sxy = SQRT(3d0)*biuni*s45/SQRT( 4d0*biuni**2-s45**2)
          ELSE
            sxy = 1d0              !assign it
         END IF
       END IF
      END IF
      !   Hill coefficients as a function of strenght ratios
      h = lb/s90/2d0
      f = 1d0/s90**2 - h
      g = 1d0        - h
      l = (4d0*f + g + h )/2d0
      m = (f + 4d0*g + h )/2d0
      n = 3d0/sxy**2/2d0
      ppr(16:21) = (/ f, g, h, l, m, n /)
      IF( .NOT.r )ppr(22:27) = ppr(16:21)
    END IF

    IF( r )THEN  ! IF Lankford ratios exist
      IF( ppr(13) /= 0d0 )THEN
        r0  = ppr(13)       ! Lankford ratio along dir 1
        r45 = ppr(14)       ! Lankford ratio for shear
        r90 = ppr(15)       ! Lankford ratio along dir 2
      ELSE
        r0  = rm             ! Uniform Lankford ratios
        r45 = rm             !
        r90 = rm             !
      END IF
      !   Hill coefficients as a function of Lankford ratios
      g = 1d0/(1d0+r0)
      h = r0*g
      f = h/r90
      l = (4d0*f + g + h )/2d0
      m = (f + 4d0*g + h )/2d0
      n = (0.5d0+r45)*(f+g)
      ppr(22:27) = (/ f, g, h, l, m, n /)
      ! recompute coefficients if both Lankford and Strengths exists
      IF( s )THEN
        ! compute average yield stress
        aux = 1d0
        aux1= 1d0
        IF( s90 /= 0d0 )THEN
          ss90 = 1d0/SQRT(f+h)
          aux = aux + s90/ss90
          aux1= aux1+ 1d0
        END IF
        IF( s45 /= 0d0)THEN
          ss45 = 2d0/SQRT(f+g+2d0*n)
          aux = aux + 2d0*s45/ss45
          aux1= aux1+ 2d0
        END IF
        IF( biuni /= 0d0 )THEN
          sb   = 1d0/SQRT(f+g)
          aux = aux+biuni/sb
          aux1= aux1+ 1d0
        END IF
        aux = (aux/aux1)**2    !average strength ratio squared
        ppr(22:27) = ppr(22:27)/aux
      ELSE                     !
        ppr(16:21) = ppr(22:27)
      END IF
    END IF
    IF( assoc == 1 ) ppr(16:21) = ppr(22:27) !associative plasticity
    WRITE(lures,"('Hill Coefs ',/,'Yield Function ',/,6f10.6,/  &
              &'Potencial Function ',/,6f10.6)",ERR=9999) ppr(16:27)

  ! new options may be computed from BIUNI and SXY for Hill-79 and Hill-90
  CASE (4)   !HILL 79
    IF( ppr(16) == 0d0 ) ppr(16) = 2d0  !Hill exponent
    IF( ppr(12) /= 0d0)THEN   !average Lankford ratio
      rm = ppr(12)
    ELSE
      IF( ppr(13) /= 0d0 )THEN
        r0  = ppr(13)       ! Lankford ratio along dir 1
        r45 = ppr(14)       ! Lankford ratio for shear
        r90 = ppr(15)       ! Lankford ratio along dir 2
        rm = (r0+2d0*r45+r90)/4d0 !average
      ELSE
        rm = 1d0
      END IF
    END IF
    ppr(17) = 0.5d0/(1d0+rm)       !K1
    ppr(18) = ppr(17)*(1d0+2d0*rm) !K2
  CASE (5)   !HILL 90
    !biuni = ppr(9)      ! biaxial strength ratio
    IF( ppr(16) == 0d0 ) ppr(16) = 2d0
    mh = ppr(16)         ! Hill exponent
    IF( ppr(13) /= 0d0 )THEN
      r0  = ppr(13)       ! Lankford ratio along dir 1
      r45 = ppr(14)       ! Lankford ratio for shear
      r90 = ppr(15)       ! Lankford ratio along dir 2
    ELSE IF( ppr(12) /= 0d0 )THEN
      r0  = ppr(12)        ! Uniform Lankford ratios
      r45 = ppr(12)        !
      r90 = ppr(12)        !
    ELSE
      r0 = 1d0
      r90= 1d0
      r45= 1d0
    END IF
    ppr(17) = 1d0+2d0*r45                                 !(sigma/tau)**m
    ppr(18) = (r0-r90)*(1d0-(mh-2d0)/2d0*r45)/((r0+r90)-(mh-2)*r0*r90) !a
    ppr(19) = mh*(r45*(r0+r90)-2d0*r0*r90)/((mh-2)*r0*r90-(r0+r90))    !b
    IF( ppr(9) > 0 )THEN
      ppr(20) = 2d0*ppr(9)
    ELSE
      ppr(20) = 2d0
    END IF
    ppr(21) = mh/2d0 - 1d0
    ppr(22) = SQRT((1d0+(1d0+r0**2)/(1+r0)**2)/1.5d0)
    ppr(23) = mh - 2d0
  END SELECT

  RETURN
  9999 CALL runen2('')
  END SUBROUTINE cohill

   SUBROUTINE coef_l(r,s90,s45,lb)
   ! IF s45 and no biaxial compute Lambda parameter
   ! I do not remember what have I done here 
   IMPLICIT NONE
   REAL(kind=8), INTENT(IN) :: r,s90,s45
   REAL(kind=8), INTENT(OUT) :: lb

   REAL(kind=8) :: s902,s452,res,dres,h,aux
   REAL(kind=8), PARAMETER :: toler = 1d-6

   s902 = 1d0/s90**2
   s452 = 4d0/s45**2
   h = 0.5d0  !initial value
   DO
     aux = 2d0*h-1d0-s902
     res = h/(h-1d0)+h/(h-s902)+(aux+s452)/aux+4d0*r+1d0
     dres= 1d0/(h-1d0)+1d0/(h-s902)+2d0/aux- &
           h/(h-1d0)**2-h/(h-s902)**2-2d0*(aux+s452)/aux**2
     h = h - res/dres
     IF( ABS(res) < toler )EXIT
   END DO
   lb = 2d0*h*s90
   RETURN
   END SUBROUTINE coef_l
