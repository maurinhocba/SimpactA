 SUBROUTINE tres10(ngaus,nnode,nelem,ndime,ntype,matno,lnods,shape,   &
                   cartd,dvolu,tresi)
 !     calculates Thermal residual
 !     for three-dimensional geometry
 USE npo_db, ONLY: iftmp,tempe
 IMPLICIT NONE

  !--- Dummy variables
  INTEGER(kind=4):: ngaus,nnode,nelem,ndime,ntype,matno(:),lnods(:,:)
  REAL(kind=8)::  shape(:,:),cartd(:,:,:,:),dvolu(:,:),tresi(:)
  !--- Local variables
  INTEGER(kind=4):: i,g,iel,isec,osec,nk,n
  REAL(kind=8):: condu,thick,f,tnod(nnode),sk(ndime),hf(nnode),dv
  TYPE(section),POINTER:: sec
  TYPE (curve), POINTER :: ck


  osec = -1          !initializes
  DO iel=1,nelem     !for each element

    isec = matno(iel)        !associated section
    IF(isec /= osec)THEN                     !if new material, get properties
      sec => psecs(isec)%p
      osec = isec                    !keep present material
      thick = 1d0
      IF( ndime == 2 .AND. ntype == 1 ) thick = sec%rprop(1)           !original thickness
      condu = sec%mtbas%prope(8)     !conductivity
      nk    = INT(sec%mtbas%props(8))        !associated curve for Heat conductivity
      CALL cur_point( nk,sec%mtbas%chead,ck) !point to the curve
    END IF

    !     zero element conductivity matrix sk(nnode,nnode)

    sk = 0d0
    hf = 0d0

    tnod (1:nnode) = tempe (1,lnods(1:nnode,iel))  !nodal temperatures

    !  Loop over all integration points

    DO g=1,ngaus
      IF( nk > 0 ) THEN
        !  Compute temperature and conductivity at integration point
        f = DOT_PRODUCT(tnod(1:nnode),shape(1:nnode,g))
        condu = inte_cr(ck%val,nk,f,i)
      END IF
      dv = dvolu(g,iel)*thick*condu
      !
      sk = MATMUL(tnod,cartd(:,:,g,iel))*dv
      hf = hf + MATMUL(cartd(:,:,g,iel),sk)

    END DO

    DO i=1,nnode                          !for each element node
      n=lnods(i,iel)                      !node number
      IF (n /= 0) THEN                    ! if node exists
        g = iftmp(1,n)                    !DOF
        IF( g > 0 )THEN                   ! if DOF exists
          tresi(g) = tresi(g) + hf(i)
        END IF
      END IF
    END DO
  END DO

RETURN
END SUBROUTINE tres10
