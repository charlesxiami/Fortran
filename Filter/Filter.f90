program main
    implicit none
    integer,parameter::n=57,t1=8,t2=10
    real,parameter::pi=3.1415926
    real x(n),y(n),h(n)
    real a,b1,b2,dt,w1,w2
    integer i
    
    open(11,file='F:\Data\Processed\ecmwf\era_interim.sat.1958-2014DJF_PC1.txt',form='formatted')
    do i=1,n
    read(11,*)x(i)
    end do
    
    dt=1.0
    w1=2*pi/10
    w2=2*pi/30
    call responseBF2(n,dt,w1,w2,a,b1,b2,h)
    call Bfilter2(n,x,y,a,b1,b2)
    
    open(21,file='F:\Data\Processed\ecmwf\era_interim.filter_10aX30a',form='binary')
    do i=1,n
    write(21)y(i)
    print*,y(i)
    end do
    
    end program main
    include 'F:\Programming\fortran\subroutine\SButterworh.f'
    