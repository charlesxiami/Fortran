	parameter(nm=56)

	real y1(nm),y2(nm)
	real y1_hf(nm),y2_hf(nm),fx(nm)

	open (14, file='PC_comb.dat')
	open (15, file='PC_comb_hf.dat')


	do kk=1,nm
	read(14,*) nn,y1(kk),y2(kk)
	end do


	  call hfilter(y1,fx,nm)

	  do kk=1,nm
	    y1_hf(kk)=fx(kk)
	  enddo

!!

	  call hfilter(y2,fx,nm)

	  do kk=1,nm
	    y2_hf(kk)=fx(kk)
	  enddo

!!



	  do kk=1,nm
	    write(15,*)kk+1956,y1_hf(kk),y2_hf(kk)
	  enddo

	



	stop
	end
	




!CCcccccccccccccccccccccccccccccccccccccccccccccccc
!CC  Fourier Transform and reverse Transform to
!CC  retain the information with 7<wavenumber<28 
!CC
        subroutine hfilter(x,fx,nm)
        dimension WKP(nm),WK(2,0:nm),x(nm),fx(nm)


        do 12 i=1,nm
12      wkp(i)=x(i)

        PI=3.1415926
        TOTAL=float(nm)

        DO 103 K=7,28
          S1=0.0
          S2=0.0
        DO 104 I=1,nm
        S1=S1+WKP(I)*COS(FLOAT(K)*FLOAT(I-1)*PI*2.0/real(nm))
104     S2=S2+WKP(I)*SIN(FLOAT(K)*FLOAT(I-1)*PI*2.0/real(nm))

        WK(1,K)=2.*S1/TOTAL
        if(k.eq.0)WK(1,K)=S1/TOTAL
        WK(2,K)=2.*S2/TOTAL
103     continue

!c   Reverse transform
        do 31 i=1,nm
        ss=0.
        do 33 k=7,28
          ss=ss+WK(1,k)*COS(FLOAT(K)*FLOAT(I-1)*PI*2.0/TOTAL)+WK(2,k)*SIN(FLOAT(K)*FLOAT(I-1)*PI*2.0/TOTAL)
33      continue
        fx(i)=ss
31      continue
        return
        end


