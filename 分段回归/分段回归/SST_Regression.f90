    program regre
    implicit none
    integer,parameter::x=180,y=45,year=38,ty=45,n1=16,n2=21,l=0
    real mam(x,y,year),ts(ty),reg1(x,y),reg2(x,y),cor1(x,y),cor2(x,y),d1(n1),d2(n2)
    real mam1(x,y,n1),mam2(x,y,n2),ts1(n1),ts2(n2)
    !ts=time series,reg=regression
    integer i,j,k
    real r,c
    
    open(10,file='F:\Data\Processed\ncep\sat\PC1_hf_sat_MAM1970-2014.dat',form='binary')
    read(10)ts
    print*,ts
    close(10)
    
    open(20,file='F:\Data\Processed\snowcover\sc.180x45.MAM1972-2009.bin',form='binary') 
    do k=1,year
    read(20)((mam(i,j,k),i=1,x),j=1,y)
    end do
    close(20)

    do j=1,y
        do i=1,x
            do k=1,n1
                mam1(i,j,k)=mam(i,j,k+l)
            end do
            do k=1,n2
                mam2(i,j,k)=mam(i,j,k+17+l)
            end do
        end do
    end do
    
    do k=1,n1
        ts1(k)=ts(k+2)
    end do
    do k=1,n2
        ts2(k)=ts(k+19)
    end do
    
    do j=1,y
        do i=1,x
            do k=1,n1
            d1(k)=mam1(i,j,k)
            end do
            r=0.0
            c=0.0
            call regression(ts1,d1,n1,r)
            call correlation(ts1,d1,n1,c)
            if ((abs(r)<5.0e+1)) then    !(abs(r)>1.0e-7).and.
            reg1(i,j)=r
            else 
            reg1(i,j)=-9.99e8
            end if 
            if ((abs(c)<1.0e+1)) then    !
                cor1(i,j)=c
            else
                cor1(i,j)=-9.99e8
            end if
        end do
    end do
    
        do j=1,y
        do i=1,x
            do k=1,n2
            d2(k)=mam2(i,j,k)
            end do
            r=0.0
            c=0.0
            call regression(ts2,d2,n2,r)
            call correlation(ts2,d2,n2,c)
            if ((abs(r)<1.0e+1)) then
            reg2(i,j)=r
            else 
            reg2(i,j)=-9.99e8
            end if 
            if ((abs(c)<1.0e+1)) then 
 !               if(c<-0.6.or.c>-0.5) then   !(abs(c)>1.0e-5).and.
                cor2(i,j)=c
            else
            cor2(i,j)=-9.99e8
  !          end if
            end if
        end do
    end do


         
    open(30,file='F:\Data\Reprocessed\mam\var_t1_PC1_hf_reg.dat',form='binary')
    write(30)reg1,cor1
    close(30)
    
    open(31,file='F:\Data\Reprocessed\mam\var_t2_PC1_hf_reg.dat',form='binary')
    write(31)reg2,cor2
    close(31)  
       
    
    end program regre
    
    include 'F:\Programming\fortran\subroutine\regression.f90'
    include 'F:\Programming\fortran\subroutine\correlation.f90'
    
    
    
    