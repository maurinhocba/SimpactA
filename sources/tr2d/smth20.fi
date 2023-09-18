 SUBROUTINE smth20(npoin,strnd,heade,eulrf)
 ! Smooths all the Gaussian variables to nodes
 IMPLICIT NONE

   !--- Dummy variables
   !    NTYPE not used, this routine only works for NTYPE =2/3 not for plane stress
   INTEGER(kind=4),INTENT(IN):: npoin
   REAL(kind=8),POINTER:: strnd(:,:) ! INTENT(OUT)
   LOGICAL, INTENT(IN) :: eulrf
   TYPE(ele20),POINTER:: heade

   !--- Constants
   REAL(kind=8),PARAMETER:: twopi=6.283185307179586476925286766559d0
   REAL(kind=8),PARAMETER:: d32  =1.5d0
   !--- Local variables
   LOGICAL        ::  plast,visco
   INTEGER(kind=4):: i, isec, osec, mtype, ipoin, &
                     nvar, nstint, ngausv, nelem, fyield
   REAL(kind=8),ALLOCATABLE:: accpn(:),      & !'ACC'umulated PoiNt volumes
                              strsg(:),      &
                              strsn(:,:)
   REAL(kind=8)            :: press,         & !volumetric stress
                              delta
   TYPE(ele20),POINTER:: elm
   TYPE(section),POINTER:: sec     !auxiliar pointer

   !Initialization
   nstint = SIZE(heade%stint)  !size of array of stresses

   nvar = 7                    !Lagrangian Formulation
   IF(eulrf) nvar = 6          !Spatial Formulation
   ALLOCATE(accpn(npoin), strsg(nvar), strsn(nvar,3), strnd(nvar,npoin))
   !initializes
   accpn(1:npoin) = 0d0            !nodal factors (inverse of element areas)
   strnd(1:nvar,1:npoin) = 0d0     !nodal internal variables
   strsg(1:nvar)  = 0d0            !local arrays
   strsn(1:nvar,1:3) = 0d0         !

   osec = -1
   !*** Bucle de elementos
   nelem = 0
   elm => heade                      !point to first element in the list
   DO                                  !loop over all the elements
     IF (.NOT.ASSOCIATED(elm)) EXIT    !last element processed => EXIT
     !check de material for set variables to transfer
     nelem = nelem + 1
     isec = elm%matno                  !material number
     IF( isec /= osec )THEN     !checks
       sec => psecs(isec)%p
       mtype  = sec%mtbas%mtype          !check material type for element
       fyield = sec%mtbas%matdef(3)      !yield function type
       plast  = fyield > 1               !plasticity
       visco = sec%mtbas%matdef(6) == 3  !Viscoplasticity
       osec  = isec
       ngausv = 0
       IF(ASSOCIATED (elm%gausv) )ngausv = SIZE(elm%gausv)  !size of array of internal variables
     END IF

    !Elemental contribution from gaussian variables to nodal variables
     IF ( mtype == 1 .AND. fyield == 2) THEN
       IF(eulrf)THEN
         delta    =(elm%gausv(1)+elm%gausv(2)+elm%gausv(4))/3d0
         strsg(1) = elm%gausv(1)- delta    !Deviatoric Finger Tensor
         strsg(2) = elm%gausv(2)- delta    !
         strsg(3) = elm%gausv(3)           !
         strsg(4) = elm%gausv(4)- delta    !
         strsg(5) = delta                  !Volum Part Finger Tensor
         strsg(6) = elm%gausv(6)           !equiv plast strain
         IF( visco ) strsg(7) = elm%gausv(5)      !equiv plast strain
       ELSE
         press    = (elm%stint(1)+elm%stint(2)+elm%stint(4))/3d0
         strsg(1) = elm%stint(1)- press    !Deviatoric Kirchhoff stress
         strsg(2) = elm%stint(2)- press    !
         strsg(3) = elm%stint(3)           !
         strsg(4) = elm%stint(4)- press    !
         strsg(5) = press                  !Volum. Mean Kirchhoff stress
         strsg(6) = SQRT(d32*(strsg(1)**2+strsg(2)**2+strsg(4)**2)+ &
                         2d0*strsg(3)**2)  !vMises
         strsg(7) = elm%gausv(6)           !equiv plast strain
       END IF
     ELSEIF(plast)THEN  !other materials and models
         press    = (elm%stint(1)+elm%stint(2)+elm%stint(4))/3d0
         strsg(1) = elm%stint(1)- press    !Deviatoric Kirchhoff stress
         strsg(2) = elm%stint(2)- press    !
         strsg(3) = elm%stint(3)           !
         strsg(4) = elm%stint(4)- press    !
         strsg(5) = press                  !Volum. Mean Kirchhoff stress
         strsg(6) = SQRT(d32*(strsg(1)**2+strsg(2)**2+strsg(4)**2)+ &
                         2d0*strsg(3)**2)  !vMises
         strsg(7) = elm%gausv(ngausv)      !equiv plast strain
     END IF

     ! extrapolation from Gauss point to nodes
     CALL smog20(nvar,elm%Area1,strsg,strsn)

     !Assemble global smoothing matrix & rhs for stresses
     DO i=1,3
       ipoin = elm%lnods(i)
       accpn(ipoin) = accpn(ipoin) + 1d0/elm%Area1
       strnd(1:nvar,ipoin) = strnd(1:nvar,ipoin) + strsn(1:nvar,i)/elm%Area1
     END DO

     elm => elm%next
   END DO

   !*** Average nodal values
   DO i=1,npoin
     IF (accpn(i) /= 0d0) THEN
       strnd(1:nvar,i) = strnd(1:nvar,i)/accpn(i)
     ELSE
       strnd(1:nvar,i) = 0d0
     END IF
   END DO

   DEALLOCATE(accpn,strsg,strsn)
 RETURN
 END SUBROUTINE smth20
