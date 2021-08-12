program climatic

implicit none
integer,parameter::year=57,lon=144,lat=73,mon=year*12+8,x=25,y=29,nm=56
real matrix(lon,lat,mon),djf(lon,lat,year),ano(lon,lat,year),ave(lon,lat),nor(lon,lat,year)
real tt(nm),ll(year),ri(nm),std(lon,lat),sx
integer i,j,t,k
!RI=regional index

open(10,file='F:\Data\Processed\ncep\sat\ncep.sat.partial.1958Jan-2015Aug',form='binary')			   
     read(10)(((matrix(i,j,t),i=1,lon),j=1,lat),t=1,mon)
close(10)
ave(lon,lat)=0

do t=1,year   
  do j=1,lat
     do i=1,lon
 !       if (abs(matrix(i,j,t))<1.0e+7) then
        djf(i,j,t)=(matrix(i,j,t*12)+matrix(i,j,t*12+1)+matrix(i,j,t*12+2))/real(3)
 !        else
!        djf(i,j,t)=djf(i,j,t)/0.0
!        end if 
!           ave(i,j)=ave(i,j)+djf(i,j,t)/year
	  end do
   end do
end do


!do j=1,lat
!    do i=1,lon
!        do k=1,year
!            ll(k)=djf(i,j,k)
!        end do
!        call getstd(ll,nm,sx)
!        std(i,j)=sx
!    end do
!end do
!
!do t=1,year
!    do j=1,lat
!        do i=1,lon
!            nor(i,j,t)=(djf(i,j,t)-ave(i,j))/std(i,j)
!        end do
!    end do
!end do

            
!!选取特定的区域求平均       
!do k=1,nm
!    tt(k)=0
!    do j=1,y
!        do i=1,x
!            tt(k)=tt(k)+nor(i+28,j+56,k+9)
!        end do
!    end do
!    ri(k)=tt(k)/(x*y)
!    print*,1956+k,ri(k)
!end do      


open(20,file='F:\Data\Processed\ncep\sat\ncep.sat.partial.1980-2014DJF',form='binary') ! DJF average for 67 years.
     do t=23,year
	    write(20)((djf(i,j,t),i=1,lon),j=1,lat) 
	 end do
close(20)


!open(40,file='F:\Data\Processed\sat\ri.nor2.dat',form='unformatted') 
!     do t=1,56
!	    write(40)ri(t)
!	 end do
!close(40)
!
!open(50,file='F:\Data\Processed\sat\ri.nor2.txt',form='formatted') 
!     do t=1,56
!	    write(50,*)ri(t)
!	 end do
!close(50)



!print*,((djf(i,j,1),i=33,57),j=41,69)

    end program climatic
    
    include "F:\Programming\fortran\subroutine\getstd.f"
