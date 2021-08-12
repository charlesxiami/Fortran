program climatic
implicit none
integer,parameter::year=57,lon=25,lat=29,n=14,p=3
real mam(lon,lat,year),ave(lon,lat),rano(lon,lat,p),rave(lon,lat,p)
real ll(n),tt(year),u(lon,lat,p),std(lon,lat),stdd,nn
integer i,j,t,k,x1,x2

open(10,file='F:\Data\Processed\ncep\sat\sat.25X29.MAM1958-2014',form='binary')
do t=1,year
read(10)((mam(i,j,t),i=1,lon),j=1,lat)
end do
close(10)
open(20,file='F:\Data\Processed\ncep\sat\sat.25X29.3p.ano',form='binary') 
open(30,file='F:\Data\Processed\ncep\sat\sat.25X29.3p.u-test',form='binary') 

do t=1,year   
  do j=1,lat
     do i=1,lon
        ave(i,j)=ave(i,j)+mam(i,j,t)/year
	  end do
   end do
end do

    do j=1,lat
        do i=1,lon
            do t=1,year
                tt(t)=mam(i,j,t)
            end do
            call getstd(tt,year,stdd)
            std(i,j)=stdd
        end do
    end do


do k=1,p
        do j=1,lat
            do i=1,lon
                do t=1,n
                rave(i,j,k)=rave(i,j,k)+mam(i,j,t+3+(k-1)*14)/real(n)
                end do
                rano(i,j,k)=rave(i,j,k)-ave(i,j)
                u(i,j,k)=rano(i,j,k)*sqrt(real(n))/std(i,j)
            end do
        end do
end do
     
    do k=1,p
        write(20)((rano(i,j,k),i=1,lon),j=1,lat)
        write(30)((u(i,j,k),i=1,lon),j=1,lat)
    end do
    
    close(20)
    close(30)
    end program climatic
    
    include "F:\Programming\fortran\subroutine\getstd.f"
