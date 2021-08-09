!******************************************************************
!     Subroutine for detect P & S arrivals on accelerogram
!     Created by Yih-Min Wu at Central Weather Bureau
!     Revised 09/10/1997
!
!     This program was test under 16 bit, +-2g, 50sps
!     3 components accelerogram.
!******************************************************************
      SUBROUTINE
     +   autopicking(xcc,np,dt,parrive,sarrive,p_wei,s_wei,is1,ren)
!------------------------------------------------------------------
!       xcc(8192,3) (input)  acclerogram unit in gal (VNE)
!       np          (input)  number of data points
!       dt          (input)  sampling interval       (sec)
!       parrive     (return) P arrival time          (sec)
!       sarrive     (return) S arrival time          (sec)
!       P_wei       (return) P arrival weighting     (sec)
!       S_wei       (return) S arrival weighting     (sec)
!------------------------------------------------------------------
	integer np

      REAL XCC(32768,3)

      INTEGER*1 itrigger
	integer i,ren

      REAL*8 TMP,SUM1,SUM2,SUM3
      real pwave_trigger,pwave_arrive,pwave_sta,pwave_lta
      real swave_trigger,swave_arrive,swave_sta,swave_lta
      integer(1) icf
      !!!!!!!!!!!!!!!!
      ren=0
      !!!!!!!!!!!!!!!!
      pwave_trigger = 2.85
      pwave_arrive = 1.25
      pwave_sta = 0.4
      pwave_lta = 40.0
      swave_trigger = 3.0
      swave_arrive = 1.25
      swave_sta = 0.5
      swave_lta = 3.0
      icf = 2
	xsum=0.0

	do nco=1,3
	  xsum=0.0
	  do i=1,1000
	    xsum=xsum+xcc(i,nco)
	  enddo
	  xsum=xsum/1000.0
	  do i=1,np
	    xcc(i,nco)=xcc(i,nco)-xsum
	  enddo
	enddo

	!-- 2000/04/24 added band pass filter before picking S arrival

      !------------------------------------------------------------
      !- change sampling interval to sampling rate
      !------------------------------------------------------------
      irate=nint(1.00001/dt)

!--------------------------------------------------------------
!-  Initial the arrivals & weightings
!--------------------------------------------------------------
      PARRIVE = 0.0
      SARRIVE = 0.0
      P_WEI   = 4.0
      S_WEI   = 4.0

!--------------------------------------------------------------
!-  Change P wave LTA & STA from seconds to sampling points
!--------------------------------------------------------------
      LTA  = NINT(PWAVE_LTA/DT)
      ISTA = NINT(PWAVE_STA/DT)

!--------------------------------------------------------------
!-  Using ISTA+100 points to calculate initial STA
!--------------------------------------------------------------
      TMP=0.0
      DO I=1,ISTA+100
         !-----------------------------------------------------
         !- Calculate characteric function
         !-----------------------------------------------------
         IF(ICF.EQ.2)THEN
           TMP=TMP+ XCC(I+1,1)*XCC(I+1,1) +
     +              (XCC(I+1,1)-XCC(I,1))*(XCC(I+1,1)-XCC(I,1))
         ELSE
           TMP=TMP+XCC(I,1)*XCC(I,1)
         ENDIF
      ENDDO
      X_STA=TMP/REAL(ISTA+100)

!--------------------------------------------------------------
!- The initial LTA is setted as STA*1.25
!--------------------------------------------------------------
      X_LTA=X_STA*1.25

!**************************************************************
!-  Start to detect P arrival, picking P arrival on V-component
!**************************************************************
      IPARR=100
      itrigger=0

      DO I=100,NP

         !-----------------------------------------------------
         !- Calculate characteric function
         !-----------------------------------------------------
         IF(ICF.EQ.2)THEN
           TMP=XCC(I,1)*XCC(I,1)+
     +         (XCC(I,1)-XCC(I-1,1))*(XCC(I,1)-XCC(I-1,1))
         ELSE
           TMP=XCC(I,1)*XCC(I,1)
         ENDIF

         !-----------------------------------------------------
         !- Update STA & LTA for each incoming data points
         !-----------------------------------------------------
         X_STA=(X_STA*(ISTA-1)+TMP)/REAL(ISTA)
         X_LTA=(X_LTA*( LTA-1)+TMP)/REAL( LTA)

         !-----------------------------------------------------
         !- Set upper limit for LTA to avoid false trigger
         !-----------------------------------------------------
         IF(X_LTA.LT.0.005)X_LTA=0.005

         !-----------------------------------------------------
         !- Calculate STA/LTA ratio for check P arrival trigger
         !-----------------------------------------------------
         TMP=X_STA/X_LTA

         !-----------------------------------------------------
         !- If STA/LTA ratio less than PWAVE_ARRIVE, keep this
         !- point, when trigger condition was met, this point
         !- define to P wave arrival.
         !-----------------------------------------------------
         IF(TMP.LE.PWAVE_ARRIVE)IPARR=I

         !-----------------------------------------------------
         !- If STA/LTA ratio bigger than PWAVE_TRIGGER, declare
         !- P wave trigger, and exit to this loop, and to go to
         !- next step to calculate P arrival's quality.
         !-----------------------------------------------------
	   if(i.lt.is1)cycle

         IF(TMP.GT.PWAVE_TRIGGER)then
           itrigger=1
           iptri=i
           exit
         ENDIF

      ENDDO

      !-----------------------------------------------------
      !- If itrigger=0, it mean P wave did not trigger
      !- no phase was picked on this station, return
      !-----------------------------------------------------
      if(itrigger.eq.0)goto 99

      !-----------------------------------------------------
      !- Change P arrival from data points to seconds
      !-----------------------------------------------------
      PARRIVE=IPARR*DT

      !-----------------------------------------------------
      !- If (P arrival + 1 second) bigger than total data
      !- length, it mean no enough data for calculating
      !- picking quality, define P arrival's weighting = 3,
      !- and return to main program
      !-----------------------------------------------------
      if((iparr+irate) .gt. np)then
        p_wei=3.0
        goto 99
      endif

!------------------------------------------------------------
!-    Calculate P arrival's quality, P arrival's quality
!-    is defined by the ratio of sum of 1 sec amplitude
!-    square after P arrival over sum of 1 sec amplitude
!-    aquare before P arrival.
!------------------------------------------------------------

      !-----------------------------------------------------
      !- calculating sum of 1 sec amplitude square
      !- after P arrival
      !-----------------------------------------------------
      sum1=0.0
      do I=iparr+1,iparr+irate
        sum1=sum1+XCC(I,1)*XCC(I,1)
      enddo
      sum1=sum1/real(irate)

      !-----------------------------------------------------
      !- calculating sum of 1 sec amplitude square
      !- before P arrival
      !-----------------------------------------------------
      sum2=0.0
      if((iparr-irate+1)<1)then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ren=1
      do I=1,iparr
        sum2=sum2+ XCC(I,1)*XCC(I,1)     
      enddo
      else
      do I=iparr-irate+1,iparr
        sum2=sum2+ XCC(I,1)*XCC(I,1)     
      enddo
      end if
      sum2=sum2/real(irate)

      !-----------------------------------------------------
      !- Calculating the ratio
      !- If Ratio  >  30 P arrival's weighting define 0
      !-    30 > R >  15 P arrival's weighting define 1
      !-    15 > R >   3 P arrival's weighting define 2
      !-     3 > R > 1.5 P arrival's weighting define 3
      !-   1.5 > R       P arrival's weighting define 4
      !-----------------------------------------------------
	if(sum2 .gt. 0.0001)then
        SUM1=sum1/sum2
	else
	  sum1=0.1
	endif
      if(sum1.ge.30.0)p_wei=0.0
      if(sum1.lt.30.0)p_wei=1.0
      if(sum1.lt.15.0)p_wei=2.0
      if(sum1.lt. 3.0)p_wei=3.0
      if(sum1.lt. 1.5)p_wei=4.0

!------------------------------------------------------------
!- Checking false trigger caused by transmission error!
!- Two type errors were occurred in our network,
!- spike    - Short time transmission error
!- DC drift - Long time transmission error
!-
!- In this program we uses the ratio of sum of 1 sec
!- amplitude square after P arrival pass 2 sec over the
!- sum of 1 sec amplitude aquare before P arrival, to
!- dectect SPIKE false trigger, if the ration less than 1.05
!- it may be SPIKE false trigger occurred, remove this P
!- arrival and return to main program.
!------------------------------------------------------------

      !-----------------------------------------------------
      !- If (P arrival + 3 sec) bigger than total data
      !- length, it mean no enough data for checking
      !- these two types false trigger. return
      !-----------------------------------------------------
      if((iparr+irate*3) .gt. np)goto 99

      !-----------------------------------------------------
      !- calculating sum of 1 sec amplitude square
      !- after P arrival pass 2 sec
      !-----------------------------------------------------
      sum3=0.0
      do I=iptri+irate*2+1,iptri+irate*3
        sum3=sum3+(XCC(I,1)*XCC(I,1))
      enddo
      sum3=sum3/real(irate)

      !-----------------------------------------------------
      !- Calculating the ratio of sum of 1 sec amplitude
      !- square after P arrival pass 2 sec over the
      !- sum of 1 sec amplitude aquare before P arrival
      !-----------------------------------------------------
	if(sum2.gt.0.001)then
        SUM3=sum3/sum2
	else
	  sum3=0.1
	endif
      !-----------------------------------------------------
      !- If ratio less than 1.05, it may be spike type
      !- false trigger occurred, remove this P arrival
      !- and return to main program.
      !-----------------------------------------------------
      IF( sum3 .lt. 1.5 )then
        parrive=0.0
        P_WEI=4.0
        S_WEI=4.0
        goto 99
      endif

!-----------------------------------------------------------
!- Detect DC drift type false trigger
!- in this program we check the 1 sec average after
!- P arrival minus 1 sec average before 1 sec of P
!- arrival If the difference bigger than 1.0 gal, it may be
!- DC drift false trigger, remove this pick and
!- return to main program
!-----------------------------------------------------------

      !-----------------------------------------------------
      !- Calculating mean after P arrival
      !-----------------------------------------------------
      sum1=0.0
      do I=iparr+1,iparr+irate
        sum1=sum1+xcc(I,1)
      enddo
      sum1=sum1/real(irate)

      !-----------------------------------------------------
      !- Calculating mean 1 sec before P arrival
      !-----------------------------------------------------
      sum2=0.0
      if((iparr-2*irate+1)<1)then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ren=1
      do I=1,iparr
        sum2=sum2+xcc(I,1)   
      enddo
      else
      do I=iparr-2*irate+1,iparr-irate
        sum2=sum2+xcc(I,1)
      enddo
      end if
      sum2=sum2/real(irate)

      !-------------------------------------------------------
      !- If the difference bigger than 1.0, it may be the DC
      !- drift false trigger occurred, remove this P arrival
      !- and return to main program.
      !-------------------------------------------------------
      IF(abs(sum1-sum2).gt.1.0)then
        parrive=0.0
        P_WEI=4.0
        S_WEI=4.0
        goto 99
      endif
!----------------***   end detect P arrival ***---------------


!**************************************************************
!-  Start to detect S arrival, picking S arrival on 2 H-ricnent
!**************************************************************

      !--------------------------------------------------------
      !- Setup picking range, picking S arrival from 2 to 42
      !- seconds after P trigger
      !--------------------------------------------------------
      ISTART= IPTRI  +  2.0/DT
      IEND  = ISTART + 40.0/DT
      IF(IEND.GT.NP)IEND=NP


      !--------------------------------------------------------
      !- Transfer S wave STA & LTA from second to points
      !--------------------------------------------------------
      LTA = NINT(SWAVE_LTA/DT)
      ISTA= NINT(SWAVE_STA/DT)

      !--------------------------------------------------------
      !- If Istart+2*ISTA bigger than total number, declare no
      !- enough signal to pick S arrival,
      !--------------------------------------------------------
      IF((ISTART+2*ista).GT.NP)goto 99


!--------------------------------------------------------------
!-  Using ISTA points to calculate initial STA
!--------------------------------------------------------------
      sum3=0.0
      DO I=ISTART+1,ISTART+ISTA

         !----------------------------------------------------
         !- using 2- horizontial component to calculate
         !- S wave characteristic function
         !----------------------------------------------------
         IF(ICF.EQ.2)THEN
           SUM1=XCC(I,2)*XCC(I,2)+(XCC(I,2)-XCC(I-1,2))*
     +                            (XCC(I,2)-XCC(I-1,2))
           SUM2=XCC(I,3)*XCC(I,3)+(XCC(I,3)-XCC(I-1,3))*
     +                            (XCC(I,3)-XCC(I-1,3))
           TMP=SUM1+SUM2
         ELSE
           TMP=XCC(I,2)*XCC(I,2)*XCC(I,2)*XCC(I,2)
     +        +XCC(I,3)*XCC(I,3)*XCC(I,3)*XCC(I,3)
         ENDIF

         sum3=sum3+TMP
      ENDDO
      X_STA=sum3/real(ista)

      !-------------------------------------------------------
      !- Initialize the LTA, is setted  = X_STA
      !-------------------------------------------------------
      X_LTA=X_STA

!-------------------------------------------------------
!- Start to Pick S wave arrival
!-------------------------------------------------------
      itrigger=0
      isarr=0

      DO I=ISTART+ISTA+1,IEND

         !---------------------------------------------------
         !- Calculating characteristic function
         !---------------------------------------------------
         IF(ICF.EQ.2)THEN
           SUM1=XCC(I,2)*XCC(I,2)+(XCC(I,2)-XCC(I-1,2))*
     +                            (XCC(I,2)-XCC(I-1,2))
           SUM2=XCC(I,3)*XCC(I,3)+(XCC(I,3)-XCC(I-1,3))*
     +                            (XCC(I,3)-XCC(I-1,3))
           TMP=SUM1+SUM2
         ELSE
           TMP=XCC(I,2)*XCC(I,2)*XCC(I,2)*XCC(I,2)
     +        +XCC(I,3)*XCC(I,3)*XCC(I,3)*XCC(I,3)
         ENDIF

         !-----------------------------------------------------
         !- Update STA & LTA for each incoming data points
         !-----------------------------------------------------
         X_STA=(X_STA*(ISTA-1)+TMP)/REAL(ISTA)
         X_LTA=(X_LTA*( LTA-1)+TMP)/REAL( LTA)

         !-----------------------------------------------------
         !- Set upper limit for LTA to avoid false trigger
         !-----------------------------------------------------
         IF(X_LTA.LT.0.05)X_LTA=0.05

         !-----------------------------------------------------
         !- Calculate STA/LTA ratio for check S arrival trigger
         !-----------------------------------------------------
         TMP=X_STA/X_LTA

         !-----------------------------------------------------
         !- If STA/LTA ratio less than SWAVE_ARRIVE, keep this
         !- point, when trigger condition was met, this point
         !- define to S wave arrival.
         !-----------------------------------------------------
         IF(TMP.LE.SWAVE_ARRIVE)ISARR=I

         !-----------------------------------------------------
         !- If STA/LTA ratio bigger than PWAVE_TRIGGER, declare
         !- P wave trigger, and exit to this loop, and to go to
         !- next step to calculate P arrival's quality.
         !-----------------------------------------------------
         IF(TMP.GT.SWAVE_TRIGGER)then
           itrigger=1
           istri=i
           exit
         ENDIF

      ENDDO

      !-----------------------------------------------------
      !- If itrigger=0, it mean S wave did not trigger
      !- no S phase picks on this station, return
      !-----------------------------------------------------
      if(itrigger.eq.0)goto 99

!-----------------------------------------------------------
!- Found S arrival
!-----------------------------------------------------------

      !-----------------------------------------------------
      !- if Isarr not settle when start picking
      !- set Isarr at 0.1 seconds before S trigger
      !-----------------------------------------------------
      IF(ISARR.EQ.0)ISARR=ISTRI-(0.1/DT)

      !-----------------------------------------------------
      !- Transfer picking result from points to seconds
      !-----------------------------------------------------
      SARRIVE=ISARR*DT

      !--------------------------------------------------------
      !- If (Istri+ 1 sec) bigger than total number, declare no
      !- enough signal to calculate S quality, set S_wei=3
      !- and return to main program
      !--------------------------------------------------------
      IF((ISTRI+irate).GT.NP)THEN
          S_WEI=3.0
          goto 99
      ENDIF

!------------------------------------------------------------
!-    Calculate S arrival's quality, S arrival's quality
!-    is defined by the ratio of sum of 1 sec amplitude
!-    square after S arrival over sum of 1 sec amplitude
!-    aquare before S arrival, using 2-horizontial component
!------------------------------------------------------------

      !-----------------------------------------------------
      !- calculating sum of 1 sec amplitude square
      !- after S arrival, using 2 horizontial component
      !-----------------------------------------------------
      sum1=0.0
      do I=isarr+1,isarr+irate
        sum1=sum1+XCC(I,2)*XCC(I,2)+XCC(I,3)*XCC(I,3)
      enddo
      SUM1=SUM1/real(irate)

      !-----------------------------------------------------
      !- calculating sum of 1 s amplitude square
      !- before S arrival, using 2 horizontial component
      !-----------------------------------------------------
      sum2=0.0
      do I=isarr-irate+1,isarr
        sum2=sum2+XCC(I,2)*XCC(I,2)+XCC(I,3)*XCC(I,3)
      enddo
      SUM2=SUM2/real(irate)


      !-----------------------------------------------------
      !- Calculate the ratio
      !- If Ratio  >  30 S arrival's weighting define 0
      !-    30 > R >  15 S arrival's weighting define 1
      !-    15 > R >   5 S arrival's weighting define 2
      !-     5 > R >   2 S arrival's weighting define 3
      !-     2 > R       S arrival's weighting define 4
      !-----------------------------------------------------
      if(sum2.gt.0.001)then
	  SUM1=sum1/sum2
	else
	  sum1=0.1
	endif
      if(sum1.ge.30.0)S_wei=0.0
      if(sum1.lt.30.0)S_wei=1.0
      if(sum1.lt.15.0)S_wei=2.0
      if(sum1.lt. 5.0)S_wei=3.0
      if(sum1.lt. 2.0)S_wei=4.0
C-------------------------------------------------------------
C     FOUND S-WAVE ARRIVE NORMAL RETURN
C-------------------------------------------------------------
99    continue
      RETURN
      END
