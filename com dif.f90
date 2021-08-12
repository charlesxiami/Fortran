program main
    implicit none
    integer,parameter::lon=144,lat=73,year=57,mon=year*12+8
    real matrix(lon,lat,mon),djf(lon,lat,year),time1(lon,lat),time2(lon,lat),time3(lon,lat)
    real diff1(lon,lat),diff2(lon,lat)
    integer i,j,t,k
    
    open(10,file='F:\Data\Process\sat_1958JAN-2015AUG',form='binary')
    open(20,file='F:\Data\Process\t1-t2',form='binary')
    open(21,file='F:\Data\Process\t3-t2',form='binary')
    read(10)(((matrix(i,j,t),i=1,lon),j=1,lat),t=1,mon)
    k=1958
    print*,k+18,k+30,k+41
    do j=1,lat 
        do i=1,lon
            do t=1,year 
                djf(i,j,t)=(matrix(i,j,t*12)+matrix(i,j,t*12+1)+matrix(i,j,t*12+2))/3.0
            end do
        end do
    end do
    do j=1,lat 
        do i=1,lon
            do t=1,12 
                time1(i,j)=time1(i,j)+djf(i,j,t+18)/12.0
            end do
        end do
    end do  
    do j=1,lat 
        do i=1,lon
            do t=1,11 
                time2(i,j)=time2(i,j)+djf(i,j,t+30)/11.0
            end do
        end do
    end do   
    do j=1,lat 
        do i=1,lon
            do t=1,14
                time3(i,j)=time3(i,j)+djf(i,j,t+41)/14.0
            end do
        end do
    end do
    do j=1,lat 
        do i=1,lon
            diff1(i,j)=time1(i,j)-time2(i,j)
            diff2(i,j)=time3(i,j)-time2(i,j)
        end do
    end do
    
    write(20)((diff1(i,j),i=1,lon),j=1,lat)
    write(21)((diff2(i,j),i=1,lon),j=1,lat)

    end program main