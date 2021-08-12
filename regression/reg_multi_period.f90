    program regre
    implicit none
    integer,parameter::x=360,y=180,year=57,n1=18,n2=26
    real mam(x,y,year),ts(year),reg1(x,y),reg2(x,y),cor1(x,y),cor2(x,y),d1(n1),d2(n2)
    real mam1(x,y,n1),mam2(x,y,n2),ts1(n1),ts2(n2)
    !ts=time series,reg=regression
    integer i,j,k
    real r,c

    open(10,file='F:\Data\Processed\ncep\sat\PC1_SAT_1958-2014MAM',form='binary')
    read(10)ts
    close(10)
    
    open(20,file='F:\Data\Processed\ncep\sst\sst.360x180.1958-2014MAM',form='binary') 
    do k=1,year
    read(20)((mam(i,j,k),i=1,x),j=1,y)
    end do
    close(20)

    do j=1,y
        do i=1,x
            do k=1,n1
                mam1(i,j,k)=mam(i,j,k+12)
            end do
            do k=1,n2
                mam2(i,j,k)=mam(i,j,k+31)
            end do
        end do
    end do
    
    do k=1,n1
        ts1(k)=ts(k+12)
    end do
    do k=1,n2
        ts2(k)=ts(k+31)
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
            if ((abs(r)>1.0e-2).and.(abs(r)<3.0e+1)) then
            reg1(i,j)=r
            else 
            reg1(i,j)=-9.99e8
            end if 
            cor1(i,j)=c
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
            if ((abs(r)>1.0e-2).and.(abs(r)<3.0e+1)) then
            reg2(i,j)=r
            else 
            reg2(i,j)=-9.99e8
            end if 
            cor2(i,j)=c
        end do
    end do

   print*,reg1
         
    open(30,file='F:\Data\Reprocessed\mam\sat_t1_PC1_reg.dat',form='binary')
    write(30)((reg1(i,j),i=1,x),j=1,y)
    close(30)
    
    open(40,file='F:\Data\Reprocessed\mam\sat_t1_PC1_cor.dat',form='binary')
    write(40)((cor1(i,j),i=1,x),j=1,y)
    close(40)
    
    open(31,file='F:\Data\Reprocessed\mam\sat_t2_PC1_reg.dat',form='binary')
    write(31)((reg2(i,j),i=1,x),j=1,y)
    close(31)
    
    open(41,file='F:\Data\Reprocessed\mam\sat_t2_PC1_cor.dat',form='binary')
    write(41)((cor2(i,j),i=1,x),j=1,y)
    close(41)
       
    
    end program regre
    
    include 'F:\Programming\fortran\subroutine\regression.f90'
    include 'F:\Programming\fortran\subroutine\correlation.f90'
    
    
    
    