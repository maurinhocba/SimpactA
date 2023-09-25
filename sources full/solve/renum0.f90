      SUBROUTINE renum0(lnods,lpntn,nelem,nnode,npoin,krenu)
!***********************************************************************
!*
!*****THIS ROUTINE SETS UP ARRAY LPNTN IF NODE RENUMBERING FOR
!*    PROFILE MINIMIZATION IS DESIRED
!*
!***********************************************************************
      IMPLICIT none
!       routine parameters
      INTEGER (kind=4),INTENT(IN) :: nelem,nnode,npoin,krenu
      INTEGER (kind=4),INTENT(IN OUT) :: lnods(:,:)
      INTEGER (kind=4),INTENT(OUT):: lpntn(:)
!       local variables
      INTEGER (kind=4) :: nstart,nlevel,nfinal,ipoin,nposi,ns,lword,    &
     &                    i,j,k,np
      INTEGER (kind=4), ALLOCATABLE :: iwork(:),iw(:,:)

      INTERFACE
        INCLUDE 'renum1.h'
        INCLUDE 'renum2.h'
        INCLUDE 'renum3.h'
      END INTERFACE

      ALLOCATE ( iw(npoin,2) )
      iw = 0
      np = 0
      DO i=1,nelem
        DO j=1,nnode
          k=lnods(j,i)
          IF(k /= 0)THEN
            IF(iw(k,1) == 0)THEN
              np = np+1
              iw(k,1)  = np     !Node number K is assigned NP
              iw(np,2) = k      !new Node NP was original Node K
            END IF
            lnods(j,i) = iw(k,1)!New node number
          END IF
        END DO
      END DO
      k = np
      DO i=1,npoin
        IF(iw(i,1) == 0)THEN    !IF node I is anassigned
          k = k+1
          iw(i,1) = k           !assign node K to unassigned I
          iw(k,2) = i           !store in K original number I
        END IF
      END DO

!.... COMPUTES THE CONNECTIONS BETWEEN DEGREES OF FREEDOM
      lword = (nnode*(nnode-1))*nelem*2 + 3*np + 1
      ALLOCATE ( iwork(lword) )
      CALL renum1(nnode,np,nelem,lnods,iwork,nposi)

!.... AUXILIAR MEMORY FOR RENUMBERING,  NODAL CONNECTIONS (NADJNT)
      nstart = 1 + nposi
      nlevel = nstart + np
      nfinal = nlevel + np - 1
      IF(nfinal > lword) THEN
        WRITE(3,"('  renumbering module requires more working space:')")
        WRITE(3,"(10x,'number of required  INTEGER (k=4) =',i10)")nfinal
        WRITE(3,"(10x,'number of allocated INTEGER (k=4) =',i10)")lword
        STOP
      END IF

!.... RENUMBERING ROUTINES
      CALL renum2(iwork(1:nstart-1),iwork(nstart:nlevel-1),             &
     &            iwork(nlevel:nfinal),lpntn,np,ns)
      CALL renum3(iwork(1:nstart-1),iwork(nstart:nlevel-1),             &
     &            iwork(nlevel:nfinal),lpntn,np,ns)

!.... OBTAIN THE INVERSE OF ARRAY LPNTN
      DO ipoin=1,npoin
        iwork(lpntn(ipoin)) = iw(ipoin,2)
      END DO

!.... PRINT RENUMBERED NODES
      IF(krenu == -1) THEN
        WRITE(3,"(//,10x,'Renumbered Nodes :',/,10x,16('-'),/,          &
     &                7x,'Old',7x,'New',17x,'New',7x,'Old',/)")
        DO ipoin=1,npoin
          lpntn(iwork(ipoin)) = ipoin
        END DO
        DO ipoin=1,npoin
          WRITE(3,"(2(5x,i5),10x,2(5x,i5))")                            &
     &             ipoin,lpntn(ipoin),ipoin,iwork(ipoin)
        END DO
      END IF

      lpntn(1:npoin) = iwork(1:npoin)

      DO i=1,nelem
        DO j=1,nnode
          k=lnods(j,i)
          IF(k /= 0) lnods(j,i) = iw(k,2) !old node number
        END DO
      END DO

      DEALLOCATE ( iwork,iw )

      RETURN

      END SUBROUTINE renum0
