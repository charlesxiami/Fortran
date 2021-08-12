C	THIS IS A PROGRAM FOR DETECTING ABRUPT CLIMATIC CHANGE 
C	BY USING MANN-KENDALL TECHNIQUE  
	PROGRAM MK 
	DIMENSION Y(1000),YY(1000),U(1000),UF(1000),UB(1000), 
     &	M(1000),MD(1000) 

        N=57
	  NYEAR=1958
C	*************************************************** 
C	* N:     SAMPLE SIZE                              * 
C	* NYEAR: FIRST YEAR OF THE TIME SERIES            * 
C	* Y(N):  ORIGINAL TIME SERIES                     * 
C	* UF(N): ORIGINAL SERIES OF U(LN)                 * 
C	* UB(N): COUNTER SERIES OF U(LN)                  * 
C	* A,B:   CRITICAL VALUE 1.96 AND -1.96            * 
C	***************************************************   	 
	OPEN(2,FILE='F:\Data\Processed\ecmwf\
     &era_interim.sat.1958-2014DJF_PC1.txt') 
	READ(2,*)(Y(I),I=1,N) 
	 do i=1,N
	  Y(i)=Y(i)/10.0
	 enddo
	CALL SMK(Y,M,MD,UF,N) 
	DO 20 I=1,N 
  20	YY(I)=Y(N+1-I) 
	CALL SMK(YY,M,MD,U,N) 
	DO 30 I=1,N 
  30	UB(I)=-U(N+1-I) 
	OPEN(3,FILE='F:\Data\Processed\ecmwf\
     &MK.ERA_interim.SAT',form='binary') 
	A=1.96 
	B=-1.96 
	DO 40 I=1,N 
	WRITE(3)UF(I),UB(I),A,B 
  50	FORMAT(1X,I4,4F8.2) 
  40	CONTINUE 
	CLOSE(3) 
	STOP 
	END 
C*********************************************************** 
	SUBROUTINE SMK(Y,M,MD,U,N) 
	DIMENSION Y(N),M(N),MD(N),U(N) 
	M(1)=0 
	DO 10 I=2,N 
	M(I)=0 
	MD(I)=0 
	DO 20 J=1,I-1 
	IF(Y(I).LT.Y(J))GOTO 20 
	M(I)=M(I)+1 
  20	CONTINUE 
	MD(I)=MD(I-1)+M(I) 
  10	CONTINUE 
	U(1)=0.0 
	DO 30 I=2,N 
	E=I*(I-1)/4.00 
	VAR=I*(I-1)*(2*I+5)/72.00 
	U(I)=(MD(I)-E)/SQRT(VAR) 
  30	CONTINUE 
	RETURN 
	END 

 