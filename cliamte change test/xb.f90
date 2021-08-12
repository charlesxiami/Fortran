!   Morlet母小波实部cos(2*3.1415926*t)*exp(-0.5*t*t)  
	program xb
    implicit none
	integer nx,nc
 	parameter (nx=63,nc=81)
!   (nc与周期的对应关系：1-2;21-4;41-8;61-16;81-32)
!   nx为时间长度
	interface
	   subroutine wavelet(input,wtr,wti,wtm,wtm2,wtm2noise,nx,nc)
	     implicit none 
	     integer nx,nc         
         real input(nx),x(-3*nx+1:4*nx)                         
         real wtr(nx,nc),wti(nx,nc),wtm(nx,nc),wtm2(nx,nc),wtm2noise(nx,nc),mywtm2(nx,nc)
	   endsubroutine wavelet
	endinterface

	integer i,j,l,r	
    real input(nx),cone(nx,nc),temp(nx,nc)
	real wtr(nx,nc),wti(nx,nc),wtm(nx,nc),wtm2(nx,nc),wtm2noise(nx,nc),mywtm2(nx,nc)


!   ----------------------------读入数据---------------------------------- 
	open(22,file='D:\Paper\study\wave\xhc.gz.txt')
	read(22,*) (input(i),i=1,nx)
    close(22) 
!   -----------------------------------------------------------------------

!   ----------------------------小波变换-----------------------------------
    call wavelet(input,wtr,wti,wtm,wtm2,wtm2noise,nx,nc)
!   -----------------------------------------------------------------------
     
!	--------------------------小波能谱输出---------------------------------
	open (22,file='D:\Paper\study\wave\xb.grd',form='binary')
	write (22) ((wtm2(i,j),j=1,nc),i=1,nx)
	close (22)

	open (22,file='D:\Paper\study\wave\myxb.grd',form='binary')
	write (22) ((mywtm2(i,j),j=1,nc),i=1,nx)
	close (22)

!   -----------------------------------------------------------------------

!   -------------------去噪声后的能谱，正值即通过显著检验------------------
    open (22,file='D:\Paper\study\wave\xb-noise.grd',form='binary')
	write (22) ((wtm2noise(i,j),j=1,nc),i=1,nx)
	close (22)
!   -----------------------------------------------------------------------

!   ----------------------------边界影响-----------------------------------
	do j=1,nc
       l=int(0.5+1.414*2/1.03*2**(float(j-1)/20.))
       r=nx-l+1
	   do i=1,l-1
	      cone(i,j)=i-l
	   enddo
	   do i=l+1,r-1
	      cone(i,j)=((i-l)*(r-i))**0.5
	   enddo
	   do i=r+1,nx
	      cone(i,j)=r-i+cone(i,j)
	   enddo
	   do i=1,nx
	      temp(i,j)=8*cone(i,j)/((nx+j)**0.5)
	   enddo	     	   	      
	enddo
	do i=2,nx-1
	do j=2,nc-1
	   cone(i,j)=(temp(i+1,j)+temp(i-1,j)+temp(i,j+1)+temp(i,j-1))*2+temp(i,j)*8
       cone(i,j)=cone(i,j)+(temp(i+1,j-1)+temp(i-1,j+1)+temp(i+1,j+1)+temp(i-1,j-1))*1
       cone(i,j)=cone(i,j)/20.
	enddo
	enddo
	do i=1,nx
	   cone(i,1)=cone(i,2)
	   cone(i,nc)=cone(i,nc-1)
	enddo
	open (22,file='D:\Paper\study\wave\cone.grd',form='binary')
	write (22)((cone(i,j),j=1,nc),i=1,nx)
	close (22)
!   ---------------------------------------------------------------------------

	end program xb




	subroutine wavelet(input,wtr,wti,wtm,wtm2,wtm2noise,nx,nc)
       implicit none 
	   integer nx,nc         
       real input(nx),x(-3*nx+1:4*nx)                         
       real wtr(nx,nc),wti(nx,nc),wtm(nx,nc),wtm2(nx,nc),wtm2noise(nx,nc),mywtm2(nx,nc)
       integer k,a,b
	   real c,noise(nc),rc
       real tba,mtr,mti,sum
	   sum=0.0                                       
  	   do b=1,nx
	      sum=sum+input(b)
       enddo
	   sum=sum/nx
	   do b=1,nx
	     input(b)=input(b)-sum
	   enddo
       sum=0.0                                      
	   do b=1,nx
	      sum=sum+input(b)**2
       enddo
	   sum=sum/nx
	   sum=sum**0.5
	   do b=1,nx
	      input(b)=input(b)/sum
	   enddo 
       sum=0.0                    
       do b=1,nx-1
	      sum=sum+input(b)*input(b+1)
	   enddo
       sum=sum/(nx-1)         	   
       rc=(1.645*((nx-2)**0.5)-1)/(nx-1) 
	   write(*,*) "lag-1自相关:",sum ," 与",rc,"比较      从而判断  红/白噪声" 
	   do b=-3*nx+1,0                     
          x(b)=0.
       enddo	 
       do b=1,nx
          x(b)=input(b)
       enddo
       do b=nx+1,4*nx
          x(b)=0.
       enddo
       do a=1,nc
	      c=2/1.03*2**(float(a-1)/20.)
	      if (abs(sum)>rc) then 
	         noise(a)=(1.0-sum*sum)/(1.0+sum*sum-2.0*sum*cos(2*3.1415926/c))
          else
	         noise(a)=1
	      endif 
          noise(a)=noise(a)*4.61/2
!      (检验 0.1－4.61   0.05－5.99    0.01－9.21)           
	      do b=1,nx
	         wtr(b,a)=0.0
	         wti(b,a)=0.0 
	         do k=b-3*c,b+3*c
	            tba=b-k
	            tba=tba/c  			 			                
                mtr=cos(6*tba)*exp(-0.5*tba*tba)/(3.14**0.25)
		        mti=sin(6*tba)*exp(-0.5*tba*tba)/(3.14**0.25)
                wtr(b,a)=wtr(b,a)+mtr*x(k)
		        wti(b,a)=wti(b,a)+mti*x(k)
   	         enddo
	         wtr(b,a)=wtr(b,a)/c**0.5
	         wti(b,a)=wti(b,a)/c**0.5
             wtm2(b,a)=wtr(b,a)**2+wti(b,a)**2
			 mywtm2(b,a)=(wtm2(b,a)*1000)/a**2         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!我用的
			 wtm(b,a)=(wtm2(b,a))**0.5		  	   
          enddo
	   enddo
	   do a=1,nc,20
	       c=2/1.03*2**(float(a-1)/20.)
	       write(*,*) c,"period noise",noise(a)
       enddo
	   do b=1,nx
       do a=1,nc
	      wtm2noise(b,a)=wtm2(b,a)-noise(a)
	   enddo   
	   enddo

   end subroutine wavelet