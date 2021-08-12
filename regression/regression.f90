    program regre
    implicit none
    integer,parameter::x=144,y=73,year=57,mon=year*12,n1=18,n2=26
    real mam(x,y,year),umatrix(x,y,year),vmatrix(x,y,year),ts(year),d(year),reg(x,y),cor(x,y)
    real mam1(x,y,n1),mam2(x,y,n2),ts1(n1),ts(n2)
    !ts=time series,reg=regression
    integer i,j,k
    real r,c

    open(10,file='F:\Data\Processed\ncep\sat\PC1_SAT_1958-2014MAM',form='binary')
    read(10)ts
    close(10)
    
    open(20,file='F:\Data\Processed\ncep\sat\sat.144X73.1958-2014MAM',form='binary') 
    do k=1,year
    read(20)((mam(i,j,k),i=1,x),j=1,y)
    end do
    close(20)

!**********************************************************************************************
    !当用到风速矢量时取消下面备注
    !
    !open(60,file='F:\Data\processed\wnd\850uwnd1957-2012.mam.dat',form='unformatted') 
    !do k=1,year
    !    read(60)((umatrix(i,j,k),i=1,x),j=1,y)
    !end do 
    !close(60)
    !open(50,file='F:\Data\processed\wnd\850vwnd1957-2012.mam.dat',form='unformatted') 
    !do k=1,year
    !    read(50)((vmatrix(i,j,k),i=1,x),j=1,y)
    !end do 
    !close(50)
    !
    ! do k=1,year
    !    do j=1,y
    !        do i=1,x
    !        mam(i,j,k)=sqrt(umatrix(i,j,k)**2+vmatrix(i,j,k)**2)
    !        end do
    !    end do
    ! end do
    ! 
!**********************************************************************************************
     
    do j=1,y
        do i=1,x
            do k=1,year
            d(k)=mam(i,j,k)
            end do
            r=0.0
            c=0.0
            call regression(ts,d,year,r)
            call correlation(ts,d,year,c)
            reg(i,j)=r
            cor(i,j)=c
        end do
    end do
!    print*,ts
         
    open(30,file='F:\Data\Reprocessed\mam\sat_PC2_reg.dat',form='binary')
    write(30)((reg(i,j),i=1,x),j=1,y)
    close(30)
    open(40,file='F:\Data\Reprocessed\mam\sat_PC2_cor.dat',form='binary')
    write(40)((cor(i,j),i=1,x),j=1,y)
    close(40)
    
    end program regre
    
    include 'F:\Programming\fortran\subroutine\regression.f90'
    include 'F:\Programming\fortran\subroutine\correlation.f90'
    
    
    
    