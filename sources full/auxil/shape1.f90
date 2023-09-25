      SUBROUTINE shape1 (nnode ,deriv ,shape ,exisp,etasp )
!********************************************************************
!
!***  calculates shape functions and their derivatives for 2d elements
!
!********************************************************************
      IMPLICIT NONE

      INTEGER (kind=4) nnode
      REAL (kind=8)    deriv(nnode,2),shape(nnode),exisp,etasp

      REAL (kind=8)    s,t,s2,t2,ss,tt,st,sst,stt,st2,stq,q2

      s  = exisp
      t  = etasp
      st = s*t


      IF(nnode == 4) THEN
        !*** shape functions for 4 noded element
        shape(1) = (1d0-s-t+st)/4d0
        shape(2) = (1d0+s-t-st)/4d0
        shape(3) = (1d0+s+t+st)/4d0
        shape(4) = (1d0-s+t-st)/4d0

        !***    and derivatives

        deriv(1,1) = (-1d0+t)/4d0
        deriv(2,1) = -deriv(1,1)
        deriv(3,1) = (1d0+t)/4d0
        deriv(4,1) = -deriv(3,1)
        deriv(1,2) = (-1d0+s)/4d0
        deriv(2,2) = (-1d0-s)/4d0
        deriv(3,2) = -deriv(2,2)
        deriv(4,2) = -deriv(1,2)

      ELSE

        s2 = s*2d0
        t2 = t*2d0
        ss = s*s
        tt = t*t
        sst = ss*t
        stt = st*t
        st2 = st*2d0


        IF (nnode == 8) THEN
          !***    shape functions for 8 noded element
          shape(1) = (-1d0+st+ss+tt-sst-stt)/4d0
          shape(2) = (-1d0-st+ss+tt-sst+stt)/4d0
          shape(3) = (-1d0+st+ss+tt+sst+stt)/4d0
          shape(4) = (-1d0-st+ss+tt+sst-stt)/4d0
          shape(5) = (1d0-t-ss+sst)/2d0
          shape(6) = (1d0+s-tt-stt)/2d0
          shape(7) = (1d0+t-ss-sst)/2d0
          shape(8) = (1d0-s-tt+stt)/2d0

          !***      and derivatives

          deriv(1,1) = (t+s2-st2-tt)/4d0
          deriv(2,1) = (-t+s2-st2+tt)/4d0
          deriv(3,1) = (t+s2+st2+tt)/4d0
          deriv(4,1) = (-t+s2+st2-tt)/4d0
          deriv(5,1) = -s+st
          deriv(6,1) = (1d0-tt)/2d0
          deriv(7,1) = -s-st
          deriv(8,1) = (-1d0+tt)/2d0
          deriv(1,2) = (s+t2-ss-st2)/4d0
          deriv(2,2) = (-s+t2-ss+st2)/4d0
          deriv(3,2) = (s+t2+ss+st2)/4d0
          deriv(4,2) = (-s+t2+ss-st2)/4d0
          deriv(5,2) = (-1d0+ss)/2d0
          deriv(6,2) = -t-st
          deriv(7,2) = (1d0-ss)/2d0
          deriv(8,2) = -t+st

        ELSE IF (nnode == 9) THEN

          ! ***     shape functions for 9 noded element

          stq = st*st
          shape(1) = (st-stt-sst+stq)/4d0
          shape(2) = (-st+stt-sst+stq)/4d0
          shape(3) = (st+stt+sst+stq)/4d0
          shape(4) = (-st-stt+sst+stq)/4d0
          shape(5) = (-t+tt+sst-stq)/2d0
          shape(6) = (s-stt+ss-stq)/2d0
          shape(7) = (t+tt-sst-stq)/2d0
          shape(8) = (-s+stt+ss-stq)/2d0
          shape(9) = 1d0-tt-ss+stq

          ! ***     and derivatives

          q2 = stt*2d0
          deriv(1,1) = (t-tt-st2+q2)/4d0
          deriv(2,1) = (-t+tt-st2+q2)/4d0
          deriv(3,1) = (t+tt+st2+q2)/4d0
          deriv(4,1) = (-t-tt+st2+q2)/4d0
          deriv(5,1) = st-stt
          deriv(6,1) = (1d0-tt+2d0*s-q2)/2d0
          deriv(7,1) = -st-stt
          deriv(8,1) = (-1d0+tt+2d0*s-q2)/2d0
          deriv(9,1) = -2.0*s+q2
          q2=sst*2.0
          deriv(1,2) = (s-st2-ss+q2)/4d0
          deriv(2,2) = (-s+st2-ss+q2)/4d0
          deriv(3,2) = (s+st2+ss+q2)/4d0
          deriv(4,2) = (-s-st2+ss+q2)/4d0
          deriv(5,2) = (-1d0+2d0*t+ss-q2)/2d0
          deriv(6,2) = -st-sst
          deriv(7,2) = (1d0+2d0*t-ss-q2)/2d0
          deriv(8,2) = st-sst
          deriv(9,2) = -2d0*t+q2

        END IF
      END IF
      RETURN
      END SUBROUTINE shape1
