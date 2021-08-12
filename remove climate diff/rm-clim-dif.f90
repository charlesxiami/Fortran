!**************************************************************************************************
!This program is to remove the climatology difference between era40 and ERA interim from ERA40 data
!Attention:由于前期数据的处理方式，ERA40的气候平均时间取35年，即1958-1992.
!**************************************************************************************************
    program remove
    implicit none
    integer,parameter::x=144,y=73,ny1=35,ny2=35,ny3=57
    !real eramon(x,y,mon1),immon(x,y,mon2),eramon2(x,y,mon1),combine_mon(x,y,mon3)
    real era(x,y,ny1),interim(x,y,ny2),era_interim(x,y,ny3)
    real era_clim(x,y),interim_clim(x,y),diff_clim(x,y)
    integer i,j,t
    
    open(11,file='F:\Data\Processed\ecmwf\era40.sat.1958jan-2002augDJF',form='binary')
    do t=1,ny1
        read(11)((era(i,j,t),i=1,x),j=1,y)
    end do
    close(11)
    open(12,file='F:\Data\Processed\ecmwf\interim.sat.1980jan-2015sepDJF',form='binary')
    do t=1,ny2
        read(12)((interim(i,j,t),i=1,x),j=1,y)
    end do
    close(12)
    
    !Calculate the climatology average of era40 and erainterim respectively
        do t=1,ny1
            do j=1,y
                do i=1,x
                    era_clim(i,j)=era_clim(i,j)+era(i,j,t)/ny1
                end  do
            end do
        end do
    
        do t=1,ny2
            do j=1,y
                do i=1,x
                    interim_clim(i,j)=era_clim(i,j)+interim(i,j,t)/ny2
                end  do
            end do
        end do

    
        do j=1,y
            do i=1,x
                diff_clim(i,j)=era_clim(i,j)-interim_clim(i,j)
            end do
        end do

    
    
    !Combine two data 
    !First period :1958jan-1979dec,ERA-40 with clim dif removed
    !Second period:1980jan-2015sep,ERA Interim
    do t=1,22
        do j=1,y
            do i=1,x
                era_interim(i,j,t)=era(i,j,t)-diff_clim(i,j)
            end do
        end do
    end do

    do t=1,ny2
        do j=1,y
            do i=1,x
                era_interim(i,j,t+22)=interim(i,j,t)
               ! print*,era_interim(1,1,t)
            end do
        end do
    end do
    
    
    !*******************************************************
    !Output combined data into two temporal resolution
    !*******************************************************
    
    open(22,file='F:\Data\Processed\ecmwf\era_interim.sat.partial.1958-2014DJF',form='binary')
    do t=1,ny3
        write(22)((era_interim(i,j,t),i=33,57),j=41,69)
    end do
    
    end program remove