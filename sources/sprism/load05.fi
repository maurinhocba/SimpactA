 SUBROUTINE load05 (igrav,force,gvect,gravy,heade,ngaus)
 !********************************************************************
 !
 !***  evaluates gravity load vector for solid element 05
 !
 !********************************************************************
 IMPLICIT NONE

 INTEGER (kind=4) :: igrav,ngaus
 REAL (kind=8) :: force(:,:),gvect(3),gravy
 TYPE (ele05), POINTER  :: heade

 !        local variable

 REAL (kind=8) :: shape(nnb,ngaus),aux(nnb,3), &
                  posgp(ngaus),weigp(ngaus),sh(3,ngaus)
 INTEGER (kind=4) i,k,n,l,isec,osec
 REAL (kind=8) weigt,fvect(3,nnb),dens,dvol
 TYPE (ele05), POINTER :: elem
 TYPE (section), POINTER :: sec
 LOGICAL :: newm
 REAL(kind=8), PARAMETER :: a = 0.3333333333333333


 IF(igrav /= 0) THEN       !if gravity is present

   !       gauss points in local coordinates and weigths
   CALL gaussq (ngaus ,posgp ,weigp )
   weigp = weigp/2d0
   k = 0
   DO i=1,ngaus
     CALL shape4(nnb,shape(1,i),aux(1,1),a,a,posgp(i))
     sh(1,i) = -0.5d0*(1d0-posgp(i))*posgp(i)
     sh(2,i) =  1d0-posgp(i)**2
     sh(3,i) =  0.5d0*(1d0+posgp(i))*posgp(i)
   END DO

   osec = -1
   elem => heade           !point to first element
   DO
     IF (.NOT.ASSOCIATED (elem) ) EXIT    !all elements processed, Exit loop

     isec  = elem%matno                     !associated material
     newm = isec /= osec
     IF ( newm )THEN
       sec => psecs(isec)%p
       IF( sec%mabas == -1 )THEN
         dens = sec%rprop(3)
       ELSE
         dens = sec%mtbas%prope(5)
       END IF
       osec = isec
       weigt = gravy*dens                     !density
     END IF
     IF(weigt /= 0.0d0) THEN                !if not zero
       fvect = 0
       DO k=1,ngaus
         dvol = DOT_PRODUCT(elem%dvol,sh(:,k))*weigp(k)
         !force vector at each node
         fvect(1,:) = fvect(1,:) +dvol*weigt*gvect(1)*shape(:,k)
         fvect(2,:) = fvect(2,:) +dvol*weigt*gvect(2)*shape(:,k)
         fvect(3,:) = fvect(3,:) +dvol*weigt*gvect(3)*shape(:,k)
       END DO
       DO n=1,nnb                                 !for each node (nnb=6)
         l = elem%lnods(n)                        !node number
         force(1:3,l) = force(1:3,l) + fvect(:,n) !sum on global force vector
       END DO
     END IF

     elem => elem%next                       !point to next element
   END DO
 END IF

 RETURN
 END SUBROUTINE load05
