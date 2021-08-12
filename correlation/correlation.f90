    program indexcorrelation
    implicit none
    integer,parameter::nm=56,year=65,mon=year*12+11
    real index1(nm),index2(nm)
    real aom(mon),aoy(year),ao(nm),yy(mon),mm(mon)
    real naom(mon),naoy(year),nao(nm),yy2(mon),mm2(mon)
    real corr1,corr2,corr3
    integer t
           
    open(10,file='F:\Data\Processed\ecmwf\sat1957-2012_PC1.dat',form='binary')
    do t=1,nm
    read(10)index1(t)
    end do
    close(10)
    open(20,file='F:\Data\Processed\sat\sat1957-2012_PC1.dat',form='unformatted')
    do t=1,nm
    read(20)index2(t) 
    end do
    close(20)
    
    open(30,file='F:\Data\Original\index\monthly.ao.index.b50.current.ascii.txt',form='formatted')
    do t=1,mon
    read(30,*)yy(t),mm(t),aom(t)
    end do
    close(30)
    open(31,file='F:\Data\Original\index\norm.nao.monthly.b5001.current.ascii.txt',form='formatted')
    do t=1,mon
    read(31,*)yy2(t),mm2(t),naom(t)
    end do
    close(31)
    
    do t=1,year
        aoy(t)=(aom(t*12)+aom(t*12+1)+aom(t*12+2))/3
        naoy(t)=(naom(t*12)+naom(t*12+1)+naom(t*12+2))/3
    end do
    
    do t=1,nm
        ao(t)=aoy(t+7)
        nao(t)=naoy(t+7)
    end do

    call correlation(index1,index2,nm,corr1)
    print*,corr1
    call correlation(index1,ao,nm,corr2)
    print*,corr2
    call correlation(index1,nao,nm,corr3)
    print*,corr3
    
    open(70,file='F:\Data\Processed\ecmwf\ec_pc1.txt')
    write(70,*)corr1
    close(70)
    !open(88,file='F:\Data\Reprocessed\work1\index_correlation\RI2.txt')
    !write(88,*)"RI2-PC1",corr1
    !write(88,*)"RI2-AO",corr2
    !write(88,*)"RI2-NAO",corr3
    !write(88,*)"RI2-cSH"
    !close(88)
    !
    !open(40,file='F:\Data\Processed\sat\ao.index.dat',form='unformatted')
    !do t=1,nm
    !write(40)ao(t)
    !end do
    !close(40)
    !open(41,file='F:\Data\Processed\sat\nao.index.dat',form='unformatted')
    !do t=1,nm
    !write(41)nao(t)
    !end do
    !close(41)
    
    end program indexcorrelation
    
    include 'F:\Programming\fortran\subroutine\correlation.f90'