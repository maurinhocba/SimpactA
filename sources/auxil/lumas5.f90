SUBROUTINE lumas5(ndofn,npoin,nelem,ngaus,nnode,lnods,matno,coord, &
     &                  emass,weigp,shape,deriv,iwrit,sumat)

! calculates lumped mass for 3d 4- or 8-node element

  USE lispa0
  USE mat_dba
  IMPLICIT NONE

  INTEGER (kind=4) ndofn,nelem,ngaus,nnode,npoin,iwrit, &
 &                 lnods(:,:),matno(:)
  REAL    (kind=8) coord(3,npoin),emass(ndofn,*), &
 &                 weigp(:),shape(:,:),deriv(:,:,:),sumat

  INTEGER (kind=4) ielem,isec,kgasp,inode,lnode,osec,istop
  REAL    (kind=8) :: djacb,dvolu,shapi,sumas,sume6,tarea,rhoel
 
  REAL(kind=8),ALLOCATABLE :: cartd(:,:,:),elcod(:,:),diagm(:)
  TYPE (section), POINTER :: sec

  ALLOCATE (cartd(nnode,3,ngaus), elcod(3,nnode), diagm(nnode))

  osec = -1
  sume6 = 0.0
  istop = 0
  DO ielem = 1,nelem
    kgasp = 0
    tarea = 0.0
    isec = matno(ielem)
    IF( isec == 0 )CYCLE     !for rigid bodies without associated node
    IF( isec /= osec )THEN
      sec => psecs(isec)%p
      ! mscal is applied over array EMASS directly
      rhoel = sec%mtbas%prope(5) !material density, (not scaled)
      osec = isec
    END IF
    IF(rhoel == 0.0d0)CYCLE        ! massless element
    DO inode=1,nnode
      diagm(1:nnode) = 0.0
      elcod(1:3,1:nnode) = coord(1:3,lnods(1:nnode,ielem))
    END DO
    DO kgasp=1,ngaus
      CALL jacob5(cartd,deriv(1,1,kgasp),djacb,elcod,nnode,ielem,istop)
      IF(istop == 1)STOP
      dvolu = djacb*weigp(kgasp)
      IF (nnode == 8) THEN
        DO inode = 1,nnode
          shapi = shape(inode,kgasp)
          diagm(inode) = diagm(inode) + shapi*dvolu !*shapi*shapi !not good for 4 gauss points
        END DO
      ELSE IF (nnode == 4) THEN ! is valid only for ngaus = 1
        DO inode = 1,nnode
          diagm(inode) = dvolu/nnode
        END DO
      ENDIF
      tarea = tarea+dvolu
    END DO

! generates lumped mass matrix proportional to diagonal

    sumas = 0.
    DO inode=1,nnode
      sumas = sumas + diagm(inode)
    END DO
    tarea = tarea*rhoel
    sume6 = sume6+tarea
    sumas = tarea/sumas
    DO inode=1,nnode
      lnode=lnods(inode,ielem)
      emass(1:3,lnode) = emass(1:3,lnode) + diagm(inode)*sumas
    END DO
  END DO
  IF(iwrit==1) WRITE(lures,"(//' Total Mass of Element 5:',e15.7)",ERR=9999) sume6
  sumat = sumat+sume6
  DEALLOCATE (cartd, elcod, diagm)

RETURN
 9999 CALL runen2('')
END SUBROUTINE lumas5
