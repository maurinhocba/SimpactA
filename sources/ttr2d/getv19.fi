 SUBROUTINE getv19(heade,elvar,eulrf)
 !-------------------------------------------------------------------
 !  Get element information to transfer
 !-------------------------------------------------------------------
 IMPLICIT NONE

   !--- Dummy arguments
   REAL(kind=8),INTENT(OUT):: elvar(:,:,:)   !Element information
   LOGICAL, INTENT(IN):: eulrf             !use finger tensor
   TYPE(ele19),POINTER:: heade
   !--- Local variables
   REAL(kind=8),PARAMETER:: d32 = 1.5d0
   REAL(kind=8) :: press,delta
   INTEGER(kind=4):: ielem, isec, osec, mtype, yield, ngausv
   TYPE(ele19),POINTER:: elm
   LOGICAL :: plast,visco
   TYPE(section),POINTER:: sec     !auxiliar pointer

   osec = -1
   ielem = 0                           !initialize element counter
   elm => heade                        !point to first element in the list
   DO                                  !loop over all the elements
     IF (.NOT.ASSOCIATED(elm)) EXIT    !last element processed => EXIT
     ielem = ielem + 1                 !increase element counter
     isec = elm%matno
     !check de material for set variables to transfer
     IF( isec /= osec )THEN     !checks
       sec => psecs(isec)%p
       mtype = sec%mtbas%mtype          !check material type for element
       yield = sec%mtbas%matdef(3)      !yield function
       plast = yield > 1                !plasticity
       visco = sec%mtbas%matdef(6) == 3 !Viscoplasticity
       osec  = isec
       IF( ASSOCIATED(elm%gausv) )THEN
         ngausv = SIZE(elm%gausv)  !size of array of internal variables
       ELSE
         ngausv = 0
       END IF
     END IF

     !Elemental contribution from gaussian variables to nodal variables
     IF ( mtype == 4 .AND. yield == 2) THEN ! yield function is MISES
       IF(eulrf)THEN
         delta    =(elm%gausv(1)+elm%gausv(2)+elm%gausv(4))/3d0
         elvar(1,1,ielem) = elm%gausv(1)- delta    !Deviatoric Finger Tensor
         elvar(2,1,ielem) = elm%gausv(2)- delta    !
         elvar(3,1,ielem) = elm%gausv(3)           !
         elvar(4,1,ielem) = elm%gausv(4)- delta    !
         elvar(5,1,ielem) = delta                  !Volum Part Finger Tensor
         elvar(6,1,ielem) = elm%gausv(6)           !equiv plast strain
         !IF( visco ) elvar(7,1,ielem) = elm%gausv(7) !plast strain rate ... not yet
       ELSE
         press = (elm%stint(1)+elm%stint(2)+elm%stint(4))/3d0
         elvar(1,1,ielem) = elm%stint(1)- press    !Deviatoric Kirchhoff stress
         elvar(2,1,ielem) = elm%stint(2)- press    !
         elvar(3,1,ielem) = elm%stint(3)           !
         elvar(4,1,ielem) = elm%stint(4)- press    !
         elvar(5,1,ielem) = press                  !Volum. Mean Kirchhoff stress
         elvar(6,1,ielem) = SQRT(d32*(elvar(1,1,ielem)**2+elvar(2,1,ielem)**2+ &
                               elvar(4,1,ielem)**2)+2d0*elvar(3,1,ielem)**2)  !vMises
         elvar(7,1,ielem) = elm%gausv(6)           !equiv plast strain
       END IF
     END IF

     elm => elm%next          !point to next node in the list
   END DO

 RETURN
 END SUBROUTINE getv19
