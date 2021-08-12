program composite_difference
implicit none
integer,parameter::year=57,lon=144,lat=73,n1=18,n2=26
real mam(lon,lat,year),time1(lon,lat),time2(lon,lat),tt(lon,lat)
real rr1(n1),rr2(n2),diff(lon,lat)
real ts
integer i,j,t,k

open(11,file='F:\Data\Processed\ncep\sat\sat.144X73.1958-2014mam',form='binary')
open(21,file='F:\Data\Processed\ncep\sat\2period_dif',form='binary') 
open(22,file='F:\Data\Processed\ncep\sat\2period_dif_Ttest',form='binary') 

do t=1,year
    read(11)((mam(i,j,t),i=1,lon),j=1,lat)
end do


    do j=1,lat 
        do i=1,lon
            do t=1,n1
                time1(i,j)=time1(i,j)+mam(i,j,t+12)/real(n1)
            end do
        end do
    end do  
    do j=1,lat 
        do i=1,lon
            do t=1,n2 
                time2(i,j)=time2(i,j)+mam(i,j,t+31)/real(n2)
            end do
        end do
    end do   

    do j=1,lat 
        do i=1,lon
            diff(i,j)=time1(i,j)-time2(i,j)
        end do
    end do

!T-test
        do j=1,lat
            do i=1,lon
                do t=1,n1
                rr1(t)=mam(i,j,t+12)
                end do
                do t=1,n2
                rr2(t)=mam(i,j,t+31)
                end do
                call ttest(rr1,rr2,n1,n2,ts)
                tt(i,j)=ts
            end do
        end do       

    write(21)((diff(i,j),i=1,lon),j=1,lat)
    write(22)((tt(i,j),i=1,lon),j=1,lat)

    
    close(11)
    close(21)
    close(22)

    end program composite_difference
    
    include "F:\Programming\fortran\subroutine\ttest.f90"

