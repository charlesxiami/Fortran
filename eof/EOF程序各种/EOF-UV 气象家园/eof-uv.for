


      PROGRAM EOF
C     THIS PROGRAM USES EOF FOR ANALYSING TIME SERIES
C     OF METEOROLOGICAL FIELD
C     M:LENTH OF TIME SERIES   !!!!!!!!!! m:时间序列长度
C     N:NUMBER OF GRID-POINTS  !!!!!!!!!! n:格点数
C     KS=-1:SELF; KS=0:DEPATURE; KS=1:STANDERDLIZED DEPATURE
C     KV:NUMBER OF EIGENVALUES WILL BE OUTPUT
C     KVT:NUMBER OF EIGENVECTORS AND TIME SERIES WILL BE OUTPUT
C     MNH=MIN(M,N)
C     EGVT=EIGENVACTORS, ECOF=TIME COEFFICIENTS FOR EGVT.
C     ER(KV,1)=LAMDA,LAMDA EIGENVALUE
C     ER(KV,2)=ACCUMULATE LAMDA
C     ER(KV,3)=THE SUM OF COMPONENTS VECTORS PROJECTED ONTO
c      EIGENVACTOR.
C     ER(KV,4)=ACCUMULATE ER(KV,3)
C

      PARAMETER(M=60,N=1218,MNH=60,KS=-1,KV=5,KVT=5,pi=3.1415926)
C
      DIMENSION F(N,M),A(MNH,MNH),S(MNH,MNH),ER(MNH,4),
     *     DF(N),V(MNH),AVF(N),EGVT(N,KVT),ECOF(M,KVT),
     *     u1(609,m),u2(609,m)
	     
	
  

      open(10,file='d:\fortran\data\u.grd',form='binary')  
      open(11,file='d:\fortran\data\v.grd',form='binary')  
	open(20,file='d:\fortran\uegvt7.grd',form='binary') 
	open(21,file='d:\fortran\vegvt7.grd',form='binary')  
	open(30,file='d:\fortran\t7.grd',form='binary') 
	open(16,file='d:\fortran\eof7.txt')    !改变数据路径

        

cccccccccccccccccc读数据
         read(10) ((u1(i,j),i=1,609),j=1,m)
	   read(11) ((u2(i,j),i=1,609),j=1,m)
cccccccccccccccc 处理数据 从hh到h的转换，再从h到f的转换
         do j=1,m
	     do i=1,609
	f(i,j)=u1(i,j)
	f(i+609,j)=u2(i,j)
	enddo;enddo


CCCCCCCCCCCCCCCCINPUT DATA CCCCCCCCCCCCCCCCCCC

      CALL TRANSF(N,M,F,AVF,DF,KS)
	write(*,*)'ok program 1'
      CALL FORMA(N,M,MNH,F,A)
	write(*,*)'ok program 2'
      CALL JCB(MNH,A,S,0.00001)
	write(*,*)'ok program 3'
      CALL ARRANG(KV,MNH,A,ER,S)
	write(*,*)'ok program 4'
      CALL TCOEFF(KVT,KV,N,M,MNH,S,F,V,ER)
	write(*,*)'ok program 5'
      CALL OUTER(KV,ER,MNH)
	write(*,*)'ok program 6'
      CALL OUTVT(KVT,N,M,MNH,S,F,EGVT,ECOF)
	write(*,*)'ok program 7'

ccccccccccccc存储数据


	do j=1,m
	do i=1,kvt
	write(30)ecof(j,i)
	enddo;enddo

      do it=1,kvt
	do j=1,609
	write(20)egvt(j,it)
	enddo;enddo	
	do it=1,kvt
	do j=610,n
	write(21)egvt(j,it)
	enddo;enddo	
      write(*,*)'ok 8'


cccccccccccc

   
      END



ccccccccccccccccccccccccc子程序

      SUBROUTINE TRANSF(N,M,F,AVF,DF,KS)
C     THIS SUBROUTINE PROVIDES INITIAL F BY KS
      DIMENSION F(N,M),AVF(N),DF(N)
      DO 5 I=1,N
      AVF(I)=0.0
  5   DF(I)=0.0
      IF(KS) 30,10,10
  10  DO 14 I=1,N
      DO 12 J=1,M
  12  AVF(I)=AVF(I)+F(I,J)
      AVF(I)=AVF(I)/M
      DO 14 J=1,M
      F(I,J)=F(I,J)-AVF(I)
  14  CONTINUE
      IF(KS.EQ.0) THEN
      RETURN
      ELSE
      DO 24 I=1,N
      DO 22 J=1,M
  22  DF(I)=DF(I)+F(I,J)*F(I,J)
      DF(I)=SQRT(DF(I)/M)
      DO 24 J=1,M
      F(I,J)=F(I,J)/DF(I)
  24  CONTINUE
      ENDIF
  30  CONTINUE
      RETURN
      END

      SUBROUTINE FORMA(N,M,MNH,F,A)
C     THIS SUBROUTINE FORMS A BY F
      DIMENSION F(N,M),A(MNH,MNH)
      IF(M-N) 40,50,50
  40  DO 44 I=1,MNH
      DO 44 J=I,MNH
      A(I,J)=0.0
      DO 42 IS=1,N
  42  A(I,J)=A(I,J)+F(IS,I)*F(IS,J)
      A(J,I)=A(I,J)
  44  CONTINUE
      RETURN
  50  DO 54 I=1,MNH
      DO 54 J=I,MNH
      A(I,J)=0.0
      DO 52 JS=1,M
  52  A(I,J)=A(I,J)+F(I,JS)*F(J,JS)
      A(J,I)=A(I,J)
  54  CONTINUE
      RETURN
      END


      SUBROUTINE JCB(N,A,S,EPS)
C     THIS SUBROUTINE COMPUTS EIGENVALUES AND standard EIGENVECTORS OF A
      DIMENSION A(N,N),S(N,N)
      DO 30 I=1,N
      DO 30 J=1,I
      IF(I-J) 20,10,20
  10  S(I,J)=1.
      GO TO 30
  20  S(I,J)=0.
      S(J,I)=0.
  30  CONTINUE
      G=0.
      DO 40 I=2,N
      I1=I-1
      DO 40 J=1,I1
  40  G=G+2.*A(I,J)*A(I,J)
      S1=SQRT(G)
      S2=EPS/FLOAT(N)*S1
      S3=S1
      L=0
  50  S3=S3/FLOAT(N)
  60  DO 130 IQ=2,N
      IQ1=IQ-1
      DO 130 IP=1,IQ1
      IF(ABS(A(IP,IQ)).LT.S3) GOTO 130
      L=1
      V1=A(IP,IP)
      V2=A(IP,IQ)
      V3=A(IQ,IQ)
      U=0.5*(V1-V3)
      IF(U.EQ.0.0) G=1.
      IF(ABS(U).GE.1E-10) G=-SIGN(1.,U)*V2/SQRT(V2*V2+U*U)
      ST=G/SQRT(2.*(1.+SQRT(1.-G*G)))
      CT=SQRT(1.-ST*ST)
      DO 110 I=1,N
      G=A(I,IP)*CT-A(I,IQ)*ST
      A(I,IQ)=A(I,IP)*ST+A(I,IQ)*CT
      A(I,IP)=G
      G=S(I,IP)*CT-S(I,IQ)*ST
      S(I,IQ)=S(I,IP)*ST+S(I,IQ)*CT
  110 S(I,IP)=G
      DO 120 I=1,N
      A(IP,I)=A(I,IP)
  120 A(IQ,I)=A(I,IQ)
      G=2.*V2*ST*CT
      A(IP,IP)=V1*CT*CT+V3*ST*ST-G
      A(IQ,IQ)=V1*ST*ST+V3*CT*CT+G
      A(IP,IQ)=(V1-V3)*ST*CT+V2*(CT*CT-ST*ST)
      A(IQ,IP)=A(IP,IQ)
  130 CONTINUE
      IF(L-1) 150,140,150
  140 L=0
      GO TO 60
  150 IF(S3.GT.S2) GOTO 50
      RETURN
      END


      SUBROUTINE ARRANG(KV,MNH,A,ER,S)
C     THIS SUBROUTINE PROVIDES A SERIES OF EIGENVALUES
C          FROM MAX TO MIN
      DIMENSION A(MNH,MNH),ER(MNH,4),S(MNH,MNH)
      TR=0.0
      DO 200 I=1,MNH
      TR=TR+A(I,I)
  200 ER(I,1)=A(I,I)
      MNH1=MNH-1
      DO 210 K1=MNH1,1,-1
      DO 210 K2=K1,MNH1
      IF(ER(K2,1).LT.ER(K2+1,1)) THEN
      C=ER(K2+1,1)
      ER(K2+1,1)=ER(K2,1)
      ER(K2,1)=C
      DO 205 I=1,MNH
      C=S(I,K2+1)
      S(I,K2+1)=S(I,K2)
      S(I,K2)=C
  205 CONTINUE
      ENDIF
  210 CONTINUE
      ER(1,2)=ER(1,1)
      DO 220 I=2,KV
      ER(I,2)=ER(I-1,2)+ER(I,1)
  220 CONTINUE
      DO 230 I=1,KV
      ER(I,3)=ER(I,1)/TR
      ER(I,4)=ER(I,2)/TR
  230 CONTINUE
      WRITE(*,250) TR
  250 FORMAT(/5X,'TOTAL SQUARE ERROR=',F20.5)
      RETURN
      END


      SUBROUTINE TCOEFF(KVT,KV,N,M,MNH,S,F,V,ER)
C     THIS SUBROUTINE PROVIDES STANDARD EIGENVECTORS (M.GE.N,SAVED IN S;
C          M.LT.N,SAVED IN F) AND ITS TIME COEFFICENTS SERIES (M.GE.N,
C          SAVED IN F; M.LT.N,SAVED IN S)
      DIMENSION S(MNH,MNH),F(N,M),V(MNH),ER(MNH,4)
	
      IF(N.LE.M) THEN
      DO 390 J=1,M
      DO 370 I=1,N
      V(I)=F(I,J)
      F(I,J)=0.
  370 CONTINUE
      DO 380 IS=1,KVT
      DO 380 I=1,N
  380 F(IS,J)=F(IS,J)+V(I)*S(I,IS)
  390 CONTINUE
      ELSE
      DO 410 I=1,N
      DO 400 J=1,M
      V(J)=F(I,J)
      F(I,J)=0.
  400 CONTINUE
      DO 410 JS=1,KVT
      DO 410 J=1,M
      F(I,JS)=F(I,JS)+V(J)*S(J,JS)
  410 CONTINUE
      DO 430 JS=1,KVT
      DO 420 J=1,M
      S(J,JS)=S(J,JS)*SQRT(ER(JS,1))
  420 CONTINUE
      DO 430 I=1,N
      F(I,JS)=F(I,JS)/SQRT(ER(JS,1))
  430 CONTINUE
      ENDIF
      RETURN
      END


      SUBROUTINE OUTER(KV,ER,MNH)
C     THIS SUBROUTINE PRINTS ARRAY ER
C     ER(KV,1) FOR  SEQUENCE OF EIGENVALUE FROM BIG TO SMALL
C     ER(KV,2) FOR  EIGENVALUE FROM BIG TO SMALL
C     ER(KV,3) FOR  SMALL LO=(LAMDA/TOTAL VARIANCE)
C     ER(KV,4) FOR  BIG LO=SUM OF SMALL LO)
      DIMENSION ER(MNH,4)
      WRITE(16,510)
  510 FORMAT(/10X,'EIGENVALUE AND ANALYSIS ERROR')
      WRITE(16,520)
  520 FORMAT(10X,1HH,8X,5HLAMDA,10X,6HSLAMDA,11X,2HPH,12X,3HSPH)
      WRITE(16,530) (IS,(ER(IS,J),J=1,4),IS=1,KV)
  530 FORMAT(1X,I10,4F15.5)
      WRITE(16,540)
  540 FORMAT(//)
      RETURN
      END


      SUBROUTINE OUTVT(KVT,N,M,MNH,S,F,EGVT,ECOF)
C     THIS SUBROUTINE PRINTS STANDARD EIGENVECTORS
C          AND ITS TIME-COEFFICENT SERIES
      DIMENSION F(N,M),S(MNH,MNH),EGVT(N,KVT),ECOF(M,KVT)
      WRITE(16,560)
  560 FORMAT(10X,'STANDARD EIGENVECTORS')
      WRITE(16,570) (IS,IS=1,KVT)
  570 FORMAT(3X,10i7)
      DO 550 I=1,N
      IF(M.GE.N) THEN
      WRITE(16,580) I,(S(I,JS),JS=1,KVT)
  580 FORMAT(1X,I3,10F7.3,/)
      DO 11 JS=1,KVT
      EGVT(I,JS)=S(I,JS)
   11 CONTINUE
      ELSE
      WRITE(16,590) I,(F(I,JS),JS=1,KVT)
  590 FORMAT(1X,I5,10F7.3)
      DO 12 JS=1,KVT
      EGVT(I,JS)=F(I,JS)
   12 CONTINUE
      ENDIF
  550 CONTINUE
C      WRITE(16,590) I,(F(I,JS),JS=1,KVT)
!      WRITE(20)((F(I,JS),i=1,n),JS=1,KVT)

      WRITE(16,720)
  720 FORMAT(//)
      WRITE(16,610)
  610 FORMAT(10X,'TIME-COEFFICENT SERIES OF S. E.')
      WRITE(16,620) (IS,IS=1,KVT)
  620 FORMAT(3X,5i12)
      DO 600 J=1,M
      IF(M.GE.N) THEN
      WRITE(16,630) J,(f(is,j),is=1,kvt)
  630 FORMAT(1X,I3,5F12.3)
      DO 13 IS=1,KVT
      ECOF(J,IS)=F(IS,J)
   13 CONTINUE
      ELSE
      WRITE(16,640) J,(S(J,IS),IS=1,KVT)
  640 FORMAT(1X,I3,10F12.3)
      DO 14 IS=1,KVT
      ECOF(J,IS)=S(J,IS)
   14 CONTINUE
      ENDIF
  600 CONTINUE
C      WRITE(30)((S(J,IS),j=1,m),IS=1,KVT)
      RETURN
      END


.