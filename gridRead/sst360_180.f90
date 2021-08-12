program climatic
implicit none
integer,parameter::year=57,lon=360,lat=180,mon=year*12
real matrix(lon,lat,mon),djf(lon,lat,year)
integer i,j,t

open(10,file='F:\Data\Original\ncep\sst.mon.1958JAN-2014DEC',form='binary')			   
     read(10)(((matrix(i,j,t),i=1,lon),j=1,lat),t=1,mon)
close(10)

  do j=1,lat
     do i=1,lon
         do t=1,year
        djf(i,j,t)=(matrix(i,j,(t-1)*12+3)+matrix(i,j,(t-1)*12+4)+matrix(i,j,(t-1)*12+5))/real(3)
        if (abs(djf(i,j,t))>1.0e+3) then 
            djf(i,j,t)=-9.99e+8
        end if
	  end do
   end do
end do

!do i=1,lon
!    do j=1,lat
!         do t=1,year    
!           ano(i,j,t)=djf(i,j,t)-ave(i,j)
!             if (abs(ano(i,j,t))>1.0e+3) then
!               ano(i,j,t)=-9.99E+8
!            end if 
!	  end do
!   end do
!end do

open(20,file='F:\Data\Processed\ncep\sst\sst.360x180.1958-2014MAM',form='binary')
     do t=1,year
	    write(20)((djf(i,j,t),i=1,lon),j=1,lat) 
	 end do
close(20)


end program climatic