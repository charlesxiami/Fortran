    program sliding_correlation
    implicit none
    integer,parameter::n=56,t=11
    real x(n),y(n),r(n-t+1)
    real xmean,ymean,sx,sy,sxy
    integer i,k
    real a,b

    !****************************
    !N:Sample size
    !T:Sliding Window
    !****************************
    open(11,file='F:\Data\Processed\sat\ri.nor2.dat',form='unformatted')
    do i=1,n
    read(11)x(i)
    end do
    close(11)
    open(12,file='F:\Data\Processed\sat\ao.index.dat',form='unformatted')
    do i=1,n
    read(12)y(i)
    end do
    close(12)
    
    do i=1,n-t+1
        xmean=0.0
        ymean=0.0
        sx=0.0
        sy=0.0
        sxy=0.0
        do k=1,t
            xmean=xmean+x(k+i-1)/t
            ymean=ymean+y(k+i-1)/t
        end do
        do k=1,t
            sxy=sxy+(x(k+i-1)-xmean)*(y(k+i-1)-ymean)
            sx=sx+(x(k+i-1)-xmean)**2
            sy=sy+(y(k+i-1)-ymean)**2
        end do
        r(i)=sxy/sqrt(sx*sy)
    end do
    
    a=0.6
    b=-0.6
    open(13,file='F:\Data\Processed\sat\ri2_ao.cor.dat',form='binary')
    do i=1,n-t+1
    write(13)r(i),a,b
    print*,1961+i,r(i)
    end do
    close(13)
    
    end program sliding_correlation