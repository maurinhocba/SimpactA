SUBROUTINE ch_d2v(flag) !,stype,curname)
!     modify curve data
USE lispa0, ONLY : lures
USE curv_db
IMPLICIT NONE

  !CHARACTER (len=*),INTENT(IN) :: stype       !not used yet
  !CHARACTER (len=mnam),INTENT(IN) :: curname  !not used yet
  INTEGER (kind=4) :: flag

  !--- Local variables
  INTEGER(kind=4)::  ltype, npts, k, ndim, l, j
  REAL(kind=8) :: var,di
  REAL(kind=8), ALLOCATABLE :: auxi(:)
  TYPE(curpar),POINTER ::  cur


  cur => tailcp !no need to seach (just added)

  ltype = INT(cur%vec(1))      !curve type (must be 3)
  IF( ltype /= 3 )THEN         !check
    WRITE(lures,"(' ERROR: prescribed displaceme must use curve defined by point')")
    WRITE(55,"(' ERROR: prescribed displaceme must use curve defined by point')")
    CALL runen2('ERROR: Prescribed displacements curve must be by points.')
  END IF
  ndim = cur%dim      !array size
  npts = (ndim-5)/2   !number of points
  ! change displacements to increments
  IF( MOD(flag,2) == 1 )THEN   !for total displacement
    l = 6              !position of first time
    di = cur%vec(l+1)  !displacement at first point
    cur%vec(l+1) = 0d0  !incremental displacement at first point
    DO k=2,npts,1       !loop over each interval
      l = l + 2       !position of present time
      var = cur%vec(l+1)  !present displacement
      cur%vec(l+1) = var - di  !incremental displacement
      di = var             !keep present displacement
    END DO
  END IF

  SELECT CASE (flag)
  CASE(1:2)  !add two points at each interval
    ALLOCATE( auxi(ndim+2) )
    auxi(1:ndim) = cur%vec
    DEALLOCATE(cur%vec)
    cur%dim = 1 + 6*npts
    ALLOCATE(cur%vec(cur%dim))
    cur%vec(1:5) = auxi(1:5)
    ! compute values at the end of the intervals
    l  = 6                                          !pointer to auxi TIME +1 => increment
    DO k=2,npts                                !for each interval
      l = l+2
      auxi(l+1) = auxi(l+1) / (auxi(l)-auxi(l-2))   !average velocity
    END DO
    l  = 6                                          !pointer to auxi TIME +1 => increment
    cur%vec(l:l+1) = auxi(l:l+1)                    !time and intial veloc
    auxi(l+1) = -auxi(l+3)
    auxi(ndim+2) = -auxi(ndim)
    j  = 8                                          !pointer to vel at the firts in-point
    DO k=2,npts                                !for each interval but the last
      l = l+2                                         !update pointer to AUXI
      var = auxi(l) - auxi(l-2)                       !interval
      cur%vec(j)   = auxi(l-2) + var*0.1d+0          !time at first point in the interval
      cur%vec(j+2) = auxi(l)   - var*0.1d+0          !time at second point in the interval
      cur%vec(j+4) = auxi(l)                          !time at the end of the interval
      var = (38d0*auxi(l+1) - auxi(l-1) - auxi(l+3))/36d0
      cur%vec(j+1) = var                              ! velocity at first point
      cur%vec(j+3) = var                              ! velocity at second point
      cur%vec(j+5) = (auxi(l+3) + auxi(l+1))/2d0      !average velocity at the end
      j = j+6                                         !update pointer to vel
    END DO
  CASE(3:4)   !add a point a each interval
    ALLOCATE( auxi(ndim) )
    auxi = cur%vec
    DEALLOCATE(cur%vec)
    ndim = 3 + 4*npts
    cur%dim = ndim
    ALLOCATE(cur%vec(ndim))
    cur%vec(1:5) = auxi(1:5)
    ! compute values at the end of the intervals
    l  = 6                                          !pointer to auxi TIME +1 => increment
    cur%vec(l) = auxi(l)                            !time
    cur%vec(l+1) = auxi(l+3)/(auxi(l+2)-auxi(l))    !average velocity
    j  = 8                                          !pointer to vel at the mid-point
    DO k=2,npts-1                                !for each interval but the last
      l = l+2                                         !update pointer to AUXI
      cur%vec(j+2) = auxi(l)                          !time at the end of the interval
      cur%vec(j) = (auxi(l-2) + auxi(l))/2d0          !time at the middle
      cur%vec(j+3) = (auxi(l+3)+auxi(l+1))/(auxi(l+2)-auxi(l-2))    !average velocity at the end
      cur%vec(j+1) = 2d0*auxi(l+1)/(auxi(l)-auxi(l-2)) - (cur%vec(j-1)+cur%vec(j+3))/2d0    !velocity at the mid-point
      j = j+4                                         !update pointer to vel
    END DO
    l = l+2                                         !update pointer to AUXI
    cur%vec(j+2) = auxi(l)                          !time at the end of the interval
    cur%vec(j) = (auxi(l-2) + auxi(l))/2d0          !time at the middle
    cur%vec(j+3) = (auxi(l+1))/(auxi(l)-auxi(l-2))    !average velocity at the end
    cur%vec(j+1) = 1.5d0*cur%vec(j+3) - cur%vec(j-1)/2d0    !velocity at the mid-point
  CASE (5:6)  !standar strategy with linear velocity variation
    ! compute velocities
    l = 6              !position of first time
    cur%vec(l+1) = 0d0 !velocity at first point
    DO k=2,npts,1       !loop over each interval
      l = l + 2       !position of present time
      var = cur%vec(l+1) /(cur%vec(l)-cur%vec(l-2))  !average velocity
      cur%vec(l+1) = 2d0*var - cur%vec(l-1) !velocity at present time
    END DO
  END SELECT
  IF(ALLOCATED(auxi) )DEALLOCATE(auxi)
RETURN
 9999 CALL runen2('')
END SUBROUTINE ch_d2v
