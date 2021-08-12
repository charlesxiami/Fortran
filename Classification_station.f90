
	implicit none
      integer,Parameter:: nm=58,ng=160,nn=2
      Real pdo(nm),nino(nm),eawm(nm),prec(ng,nm),comp(ng,nn),lev(nn),lat(ng),lon(ng),x(nm),y(nm),tt(ng,nn),tn(3),dhl
      integer i,k,j,nflag,nlev,nh,nl,mh,ml
      character*8 id(ng),num(ng)
      real time

      open(10,file='prec1957-2015DJF.dat',form='unformatted')
      open(22,file='NINO3.4index1958-2015DJF.txt')


      open(100,file="station.dat")    !读取站点资料    
             do i=1,ng
                read(100,*)num(i),id(i),lat(i),lon(i)
             enddo
	  
      read(10)(prec(k,1),k=1,ng)
      do i=1,nm
         read(10)(prec(k,i),k=1,ng)
         read(22,*)nino(i)
      end do
      print*,nino
      

      
 do k=1,ng
 !     write(30,*)'1.pdo>0
      do i=1,nm
         if((((i+1957)>=1958.and.(i+1957)<=1960).or.((i+1957)>=1978.and.(i+1957)<=2005).or.((i+1957)>=2013.and.(i+1957)<=2015)))then
            x(i)=nino(i)
            y(i)=prec(k,i)
         else
              x(i)=0.
         end if
     !    if(nino(i)>0.8.or.nino(i)<-0.8)x(i)=0.
      end do
      call differencehl1(nm,x,y,0.8,-0.8,dhl,nh,nl,tn)
      comp(k,1)=dhl                 !case1
      tt(k,1)=tn(3)            


!      write(30,*)'2.pdo<0
      do i=1,nm
         if((((i+1957)>=1961.and.(i+1957)<=1977).or.((i+1957)>=2006.and.(i+1957)<=2012)))then
            x(i)=nino(i)
            y(i)=prec(k,i)
         else
              x(i)=0.
         end if
   !      if(nino(i)>0.8.or.nino(i)<-0.8)x(i)=0.
      end do
      call differencehl1(nm,x,y,0.8,-0.8,dhl,mh,ml,tn)
      comp(k,2)=dhl                 !case2
      tt(k,2)=tn(3) 
end do
print*,nh,nl,mh,ml


  !输出结果***************************************** 
     open(40,file='E:/PDO/Prec-PDO_0.8.bin',form='binary')

      time=0.0                      
      nlev=1
      nflag=1       
      do k=1,ng
         write(40)id(k),lat(k),lon(k),time,nlev,nflag,comp(k,1),comp(k,2),tt(k,1),tt(k,2)
      end do
      nlev=0
      write(40)id(ng),lat(ng),lon(ng),time,nlev,nflag,comp(ng,1),comp(ng,2),tt(ng,1),tt(ng,2)

      end




!c-----*----------------------------------------------------6---------7--
!c     For a time series f(n), computing:
!c       (1) mean in the high index years of x(n)
!c       (2) mean in the low index years of x(n)
!c       (3) composite difference between the mean of f(n) in high 
!c                     index years of x(n) and its climatology average
!c       (4) composite difference between the mean of f(n) in low 
!c                     index years of x(n) and its climate average
!c       (5) composite difference between the means of f(n) in high 
!c                     and low index years of x(n)
!c     input: n,x(n),f(n),coefh,coefl
!c       n: the length of time series
!c       x: control variable (index)
!c       f: given series
!c       coefh: control parameter for high index (i.e., high index years are 
!c              those in which x(i) > coefh)
!c       coefl: control parameter for low index (i.e., low index years are 
!c              those in which x(i) < coefl)
!c     output: fh,fl,dh,dl,dhl,tn(5)
!c       fh: the mean of f in high index years of x(n)
!c       fl: the mean of f in low index years of x(n)
!c       dh: composite difference between the mean of f in high index years of x(n) 
!c            and its climate mean (i.e., high index years minus cliamte mean).
!c       dl: composite difference between the mean of f in low index years of x(n) 
!c            and its climate mean (i.e., low index years minus cliamte mean). 
!c       dhl: composite difference between the means of f in high and low index years
!c             of x(n) (i.e., high minus low index years) 
!c       tn(i,j): tn only equals 2., -2., 1., -1. or 0. corresponding to significant difference or not.
!c           tn=2. indicates that the difference is positive and significant.
!c           tn=-2. indicates that the difference is negative and significant.
!c           tn=1. indicates that the difference is positive but not significant.
!c           tn=-1. indicates that the difference is negative but not significant.
!c           tn=0. indicates the difference is zero.
!c           tn(1,j)~tn(5,j) are corresponding to the 90%,95%,98%,99% and 99.9% confident levels.
!c           j=1: siginificant test for dh
!c           j=2: siginificant test for dl
!c           j=3: siginificant test for dhl
!c     Feburary 11, 2002 by Jianping Li.
      subroutine differencehl1(n,x,f,coefh,coefl,dhl,nh,nl,tn)
      dimension x(n),f(n),tn(3)
      real hn(n),ln(n) !work array
      call meanvar(n,f,avef,sf,vf)
      fh=0.
      fl=0.
      nh=0
      nl=0
      do k=1,n
        if(x(k).gt.coefh)then
          nh=nh+1
          hn(nh)=f(k)
        endif
        if(x(k).lt.coefl)then
          nl=nl+1
          ln(nl)=f(k)
        endif
      enddo
      call meanvar1(nh,hn,fh)
      call meanvar1(nl,ln,fl)
      dh=fh-avef
      dl=fl-avef
      dhl=fh-fl
      call diff_t_test(nh,n,hn,f,tn(1))
      call diff_t_test(nl,n,ln,f,tn(2))
      call diff_t_test(nh,nl,hn,ln,tn(3))
      return
      end
! calculate the t test statistic values**********************
subroutine diff_t_test(n,m,x,y,ta)
      parameter(nn=10000)
      real x(n),y(m),tn
      real ft(nn,5) !work array
        tn=0.
      call meanvar(n,x,ax,sx,vx)
      call meanvar(m,y,ay,sy,vy)
      sxy=(n*vx+m*vy)/(n+m-2.)
      sxy=sxy*(1./n+1./m)
      sxy=sqrt(sxy)
      if(sxy.eq.0.)return
      dxy=ax-ay
      ta=abs(dxy)/sxy
      if(dxy.lt.0.)ta=ta*(-1)
      return

      end
!c-----*----------------------------------------------------6---------7--
!c     Computing the mean ax, standard deviation sx
!c       and variance vx of a series x(i) (i=1,...,n).
!c     input: n and x(n)
!c       n: number of raw series
!c       x(n): raw series
!c     output: ax, sx and vx
!c       ax: the mean value of x(n)
!c       sx: the standard deviation of x(n)
!c       vx: the variance of x(n)
!c     By Dr. LI Jianping, May 6, 1998.
      subroutine meanvar(n,x,ax,sx,vx)
      dimension x(n)
      ax=0.
      vx=0.
      sx=0.
      do 10 i=1,n
        ax=ax+x(i)
  10  continue
      ax=ax/float(n)
      do 20 i=1,n
        vx=vx+(x(i)-ax)**2
  20  continue
      vx=vx/float(n)
      sx=sqrt(vx)
      return
      end
!c-----*----------------------------------------------------6---------7--
!c     Computing the mean ax of a series x(i) (i=1,...,n).
!c     input: n and x(n)
!c       n: number of raw series
!c       x(n): raw series
!c     output: ax
!c       ax: the mean value of x(n)
!c     By Dr. LI Jianping, May 6, 1999.
      subroutine meanvar1(n,x,ax)
      dimension x(n)
      ax=0.
      vx=0.
      sx=0.
      do 10 i=1,n
        ax=ax+x(i)
  10  continue
      ax=ax/float(n)
      return
      end
