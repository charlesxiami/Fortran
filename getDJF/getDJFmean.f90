    program getdjf
    implicit none
    integer,parameter::lon=144,lat=73,year=44,mon=year*12+8
    real matrix(lon,lat,mon),djf(lon,lat,year),djf2(lon,lat,year)
    integer i,j,t
    open(11,file='F:\Data\Original\ecmwf\era40\era40.sat.1958jan-2002aug.dat',form='binary')
    read(11)(((matrix(i,j,t),i=1,lon),j=1,lat),t=1,mon)
    close(11)
    do t=1,year   
        do j=1,lat
            do i=1,lon
                djf(i,j,t)=(matrix(i,j,t*12)+matrix(i,j,t*12+1)+matrix(i,j,t*12+2))/real(3)
               djf2(i,j,t)=djf(i,j,t)-273.15
            end do
        end do
    end do
    open(21,file='F:\Data\Processed\ecmwf\era40.sat.partial.1958-2002DJF',form='binary')
    do t=1,year
            write(21)((djf2(i,j,t),i=33,57),j=41,69)
        end do
    close(21)
    end program getdjf
    