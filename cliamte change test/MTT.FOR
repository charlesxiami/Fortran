C	THIS IS A PROGRAM FOR DETECTING ABRUPT CLIMATIC CHANGE
C	BY USING MOVING t-TEST TECHNIQUE
	PROGRAM MTT
	DIMENSION Y(1000),YN1(1000),YN2(1000),YC1(1000),YC2(1000),
     &	S(1000),T(1000),NY1(1000),NY2(1000),NNY1(1000),NNY2(1000),
     &     x(1000)
      
	!WRITE(*,10)
 ! 10	FORMAT(5X,'N=?,IH=?,NYEAR=?')
	N=44
      IH=5
      NYEAR=1958
C	**************************************************
C	* N:    SAMPLE SIZE                              *
C	* IH:   LENGTH OF SUB-SERIES	                 *
C	* NYEAR: FIRST YEAR OF THE TIME SERIES           *
C	* Y(N):  ORIGINAL TIME SERIES                    *
C	* ************************************************
	OPEN(2,FILE='F:\Data\Processed\ecmwf\
     &PC1_era40.sat.1958-2002DJF',form='binary')
      do i=1,n
	READ(2)Y(I)
      end do
	N1=N-IH+1
	N2=N-2*IH+1
	C1=0.0
	C2=0.0
	DO 20 I=1,IH
  20	C1=C1+Y(I)
	DO 30 I=IH+1,2*IH
  30	C2=C2+Y(I)
	DO 40 I=1,N-IH
	D1=C1-Y(I)+Y(I+IH)
	YN1(I)=C1/IH
  40	C1=D1
	DO 50 I=IH+1,N-IH
	D2=C2-Y(I)+Y(I+IH)
	YN2(I-IH)=C2/IH
  50	C2=D2
	YN1(N1)=C1/IH
	YN2(N2)=C2/IH
	DO 60 I=1,N1
	YC1(I)=0.0
	DO 70 J=I,IH+I-1
  70	YC1(I)=YC1(I)+(Y(J)-YN1(I))*(Y(J)-YN1(I))
  60	CONTINUE
	DO 80 I=1,N2
	YC2(I)=0.0
	DO 90 J=IH+I,2*IH+I-1
  90	YC2(I)=YC2(I)+(Y(J)-YN2(I))*(Y(J)-YN2(I))
  80	CONTINUE
	DO 100 I=1,N2
  100	S(I)=SQRT((YC1(I)+YC2(I))/(IH+IH-2))
	DO 110 I=1,N2
  110	T(I)=(YN1(I)-YN2(I))/(S(I)*SQRT(2.0/IH))
	DO 120 I=1,N2
	NY1(I)=NYEAR+I-1
	NNY1(I)=NY1(I)+IH-1
	NY2(I)=NYEAR+I+IH-1
  120	NNY2(I)=NY2(I)+IH-1
	IF(IH.EQ.5)THEN
	A=2.31
	B=-2.31
	ELSE IF(IH.EQ.10)THEN 
	A=2.10
	B=-2.10
	ELSE
	A=3.0
	B=-3.0
	END IF
	!OPEN(3,FILE='RI2_5a.txt ',STATUS='NEW')
	!WRITE(3,130)
 ! 130	FORMAT(30X,'ABRUPT CLIMATIC CHANGE ANALYSIS'/)
	!WRITE(3,140)
 ! 140	FORMAT(25X,'TM',4X,'t-TEST(0.01)'/)
	!DO 150 I=1,N2
	!WRITE(3,160)NY1(I),NNY1(I),NY2(I),NNY2(I),T(I),A,B
 ! 160	FORMAT(1X,I4,'-',I4,'--',I4,'-',I4,1X,3F8.2)
 ! 150	CONTINUE
	OPEN(4,FILE='F:\Data\Processed\ecmwf\
     &MTT.PC1.era40.sat.1958-2002DJF.txt')
	DO 170 I=1,N2
	WRITE(4,*)NNY1(I),T(I),A,B
  180	FORMAT(1X,I4,3F8.2)
  170 CONTINUE
      open(5,file='F:\Data\Processed\ecmwf\
     &MTT.PC1.era40.sat.1958-2002DJF',form='binary')
      DO I=1,N2
	WRITE(5)T(I),A,B
      end do
      
	END
