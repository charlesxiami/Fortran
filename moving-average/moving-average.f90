
    program main
    implicit none
    integer,parameter::n=56,t=11
    real x(n),y(n),ri(n)
    integer i
    
    open(10,file='F:\Data\Processed\sat\ri.nor2.txt',form='formatted')
    do i=1,n
    read(10,*)x(i)
    end do
    
    call  runningmean(x,y,n,t)
    
    open(11,file='F:\Data\Processed\sat\ri.nor2.dat',form='unformatted') 
    do i=1,n
        read(11)ri(i)
    end do
    close(11)
    
    open(20,file='F:\Data\Processed\sat\ri2.11a.dat',form='binary')
    do i=1,n
        write(20)ri(i),y(i)
    end do
    open(21,file='F:\Data\Processed\sat\ri2.11a.txt',form='formatted')
    do i=1,n
        write(21,*)y(i)
    end do
    
    close(20)
    close(21)

    end program main
    include 'F:\Programming\fortran\subroutine\runningMEAN.f90'

