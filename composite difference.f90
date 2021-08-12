program composite_difference
implicit none
integer,parameter::year=57,lon=144,lat=73,n=12
real djf(lon,lat,year),time1(lon,lat),time2(lon,lat),time3(lon,lat),tt(lon,lat),tt2(lon,lat)
real rr1(12),rr2(12),rr3(14),diff1(lon,lat),diff2(lon,lat)
real ts,ts2
integer i,j,t

open(11,file='F:\Data\Processed\ncep\sat\sat.144X73.1958-2014DJF',form='binary')
open(21,file='F:\Data\Processed\ncep\sat\t1-t2_dif',form='binary') 
open(22,file='F:\Data\Processed\ncep\sat\t3-t2_dif',form='binary') 
open(23,file='F:\Data\Processed\ncep\sat\t1-t2_ttest',form='binary') 
open(24,file='F:\Data\Processed\ncep\sat\t3-t2_ttest',form='binary') 

do t=1,year
    read(11)((djf(i,j,t),i=1,lon),j=1,lat)
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

!T-test
        do j=1,lat
            do i=1,lon
                do t=1,12
                rr1(t)=djf(i,j,t+18)
                rr2(t)=djf(i,j,t+18+12)
                end do
                call ttest(rr1,rr2,n,n,ts)
                tt(i,j)=ts
            end do
        end do 
        
        do j=1,lat
            do i=1,lon
                do t=1,11
                rr2(t)=djf(i,j,t+30)
                rr3(t)=djf(i,j,t+41)
                end do
                call ttest(rr2,rr3,11,11,ts2)
                tt2(i,j)=ts2
            end do
        end do            

    write(21)((diff1(i,j),i=1,lon),j=1,lat)
    write(22)((diff2(i,j),i=1,lon),j=1,lat)
    write(23)((tt(i,j),i=1,lon),j=1,lat)
    write(24)((tt2(i,j),i=1,lon),j=1,lat)
    
    
    close(11)
    close(21)
    close(22)

    end program composite_difference
    
    include "F:\Programming\fortran\subroutine\ttest.f90"

