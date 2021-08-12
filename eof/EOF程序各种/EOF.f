c-----*----------------------------------------------------6---------7--
C     EMPIRICAL ORTHOGONAL FUNCTIONS (EOF's)
c     This subroutine applies the EOF approach to analysis time series 
c       of meteorological field f(m,n).
c     input: m,n,mnl,f(m,n),ks
c       m: number of grid-points
c       n: lenth of time series
c       mnl=min(m,n)
c       f(m,n): the raw spatial-temporal seires
c       ks: contral parameter
c           ks=-1: self, i.e., for the raw time series; 
c           ks=0: departure, i.e., for the anomaly time series from climatological mean;
c           ks=1: normalized departure, i.e., for the normalized time series; 
c     output: egvt,ecof,er
c       egvt(m,mnl): array of eigenvactors
c       ecof(mnl,n): array of time coefficients for the respective eigenvectors
c       er(mnl,1): lamda (eigenvalues), its sequence is from big to small.
c       er(mnl,2): accumulated eigenvalues from big to small
c       er(mnl,3): explained variances (lamda/total explain) from big to small
c       er(mnl,4): accumulated explaned variances from big to small
c     work variables:

c     Last updated by Dr. Jianping Li on October 20, 2005.

c     Citation: Li, J., 2005: EOF subroutine, http://web.lasg.ac.cn/staff/ljp/subroutine/EOF.f.

c-----
      subroutine eof(m,n,mnl,f,ks,er,egvt,ecof)
      dimension f(m,n),er(mnl,4),egvt(m,mnl),ecof(mnl,n)
      dimension cov(mnl,mnl),s(mnl,mnl),d(mnl),v(mnl) !work array
c---- Preprocessing data
      call transf(m,n,f,ks)
c---- Crossed product matrix of the data f(m,n)
      call crossproduct(m,n,mnl,f,cov)
c---- Eigenvalues and eigenvectors by the Jacobi method 
      call jacobi(mnl,cov,s,d,0.00001)
c---- Specified eigenvalues 
      call arrang(mnl,d,s,er)
c---- Normalized eignvectors and their time coefficients 
      call tcoeff(m,n,mnl,f,s,er,egvt,ecof)
      return
      end
c-----*----------------------------------------------------6---------7--
c     Preprocessing data to provide a field by ks.
c     input: m,n,f
c       m: number of grid-points
c       n: lenth of time series
c       f(m,n): the raw spatial-temporal seires
c       ks: contral parameter
c           ks=-1: self, i.e., for the raw time series; 
c           ks=0: departure, i.e., for the anomaly time series from climatological mean;
c           ks=1: normalized departure, i.e., for the normalized time series; 
c     output: f
c       f(m,n): output field based on the control parameter ks.
c     work variables: fw(n)
      subroutine transf(m,n,f,ks)
      dimension f(m,n)
      dimension fw(n),wn(m)           !work array
	i0=0
	do i=1,m
	  do j=1,n
          fw(j)=f(i,j)
        enddo
	  call meanvar(n,fw,af,sf,vf)
	  if(sf.eq.0.)then
	    i0=i0+1
	    wn(i0)=i
	  endif
	enddo
	if(i0.ne.0)then
	  write(*,*)'****  FAULT  ****'
	  write(*,*)' The program cannot go on because 
     *the original field has invalid data.'
	  write(*,*)' There are totally ',i0,
     *    '  gridpionts with invalid data.'
     	  write(*,*)' These invalid data are those whose standard variance 
     *equal zero.'
     	  write(*,*)' The array WN stores the positions of those invalid 
     *grid-points. You must pick off those invalid data from the orignal
     *field and then reinput a new field to calculate its EOFs.'   
	  write(*,*)'****  FAULT  ****'
	  stop
	endif	    
      if(ks.eq.-1)return
      if(ks.eq.0)then                !anomaly of f
        do i=1,m
          do j=1,n
            fw(j)=f(i,j)
          enddo
          call meanvar(n,fw,af,sf,vf)
          do j=1,n
            f(i,j)=f(i,j)-af
          enddo
        enddo
        return
      endif
      if(ks.eq.1)then                 !normalizing f
        do i=1,m
          do j=1,n
            fw(j)=f(i,j)
          enddo
          call meanvar(n,fw,af,sf,vf)
          do j=1,n
            f(i,j)=(f(i,j)-af)/sf
          enddo
        enddo
      endif
      return
      end
c-----*----------------------------------------------------6---------7--
c     Crossed product martix of the data. It is n times of 
c       covariance matrix of the data if ks=0 (i.e. for anomaly). 
c     input: m,n,mnl,f
c       m: number of grid-points
c       n: lenth of time series
c       mnl=min(m,n)
c       f(m,n): the raw spatial-temporal seires
c     output: cov(mnl,mnl)  
c       cov(m,n)=f*f' or f'f dependes on m and n.
c         It is a mnl*mnl real symmetric matrix.
      subroutine crossproduct(m,n,mnl,f,cov)
      dimension f(m,n),cov(mnl,mnl)
      if(n-m) 10,50,50
  10  do 20 i=1,mnl
      do 20 j=i,mnl
        cov(i,j)=0.0
        do is=1,m
          cov(i,j)=cov(i,j)+f(is,i)*f(is,j)
        enddo
        cov(j,i)=cov(i,j)
  20  continue
      return
  50  do 60 i=1,mnl
      do 60 j=i,mnl
        cov(i,j)=0.0
        do js=1,n
          cov(i,j)=cov(i,j)+f(i,js)*f(j,js)
        enddo
        cov(j,i)=cov(i,j)
  60  continue
      return
      end
c-----*----------------------------------------------------6---------7--
c     Computing eigenvalues and eigenvectors of a real symmetric matrix
c       a(m,m) by the Jacobi method.
c     input: m,a,s,d,eps 
c       m: order of matrix
c       a(m,m): the covariance matrix
c       eps: given precision
c     output: s,d 
c       s(m,m): eigenvectors
c       d(m): eigenvalues
      subroutine jacobi(m,a,s,d,eps)
      dimension a(m,m),s(m,m),d(m)
      do 30 i=1,m
      do 30 j=1,i
        if(i-j) 20,10,20
  10    s(i,j)=1.
        go to 30
  20    s(i,j)=0.
        s(j,i)=0.
  30  continue
      g=0.
      do 40 i=2,m
        i1=i-1
        do 40 j=1,i1
  40      g=g+2.*a(i,j)*a(i,j)
      s1=sqrt(g)
      s2=eps/float(m)*s1
      s3=s1
      l=0
  50  s3=s3/float(m)
  60  do 130 iq=2,m
        iq1=iq-1
        do 130 ip=1,iq1
        if(abs(a(ip,iq)).lt.s3) goto 130
        l=1
        v1=a(ip,ip)
        v2=a(ip,iq)
        v3=a(iq,iq)
        u=0.5*(v1-v3)
        if(u.eq.0.0) g=1.
        if(abs(u).ge.1e-10) g=-sign(1.,u)*v2/sqrt(v2*v2+u*u)
        st=g/sqrt(2.*(1.+sqrt(1.-g*g)))
        ct=sqrt(1.-st*st)
        do 110 i=1,m
          g=a(i,ip)*ct-a(i,iq)*st
          a(i,iq)=a(i,ip)*st+a(i,iq)*ct
          a(i,ip)=g
          g=s(i,ip)*ct-s(i,iq)*st
          s(i,iq)=s(i,ip)*st+s(i,iq)*ct
  110     s(i,ip)=g
        do 120 i=1,m
          a(ip,i)=a(i,ip)
  120     a(iq,i)=a(i,iq)
        g=2.*v2*st*ct
        a(ip,ip)=v1*ct*ct+v3*st*st-g
        a(iq,iq)=v1*st*st+v3*ct*ct+g
        a(ip,iq)=(v1-v3)*st*ct+v2*(ct*ct-st*st)
        a(iq,ip)=a(ip,iq)
  130 continue
      if(l-1) 150,140,150
  140 l=0
      go to 60
  150 if(s3.gt.s2) goto 50
      do 160 i=1,m
        d(i)=a(i,i)
  160 continue
      return
      end
c-----*----------------------------------------------------6---------7--
c     Provides a series of eigenvalues from maximuim to minimuim.
c     input: mnl,d,s
c       d(mnl): eigenvalues 
c       s(mnl,mnl): eigenvectors
c     output: er
c       er(mnl,1): lamda (eigenvalues), its equence is from big to small.
c       er(mnl,2): accumulated eigenvalues from big to small
c       er(mnl,3): explained variances (lamda/total explain) from big to small
c       er(mnl,4): accumulated explaned variances from big to small
      subroutine arrang(mnl,d,s,er)
      dimension d(mnl),s(mnl,mnl),er(mnl,4)
      tr=0.0
      do 10 i=1,mnl
        tr=tr+d(i)
        er(i,1)=d(i)
  10  continue
      mnl1=mnl-1
      do 20 k1=mnl1,1,-1
      do 20 k2=k1,mnl1
        if(er(k2,1).lt.er(k2+1,1)) then
          c=er(k2+1,1)
          er(k2+1,1)=er(k2,1)
          er(k2,1)=c
          do 15 i=1,mnl
            c=s(i,k2+1)
            s(i,k2+1)=s(i,k2)
            s(i,k2)=c
  15      continue
        endif
  20  continue
      er(1,2)=er(1,1)
      do 30 i=2,mnl
        er(i,2)=er(i-1,2)+er(i,1)
  30  continue
      do 40 i=1,mnl
        er(i,3)=er(i,1)/tr
        er(i,4)=er(i,2)/tr
  40  continue
      return
      end
c-----*----------------------------------------------------6---------7--
c     Provides standard eigenvectors and their time coefficients
c     input: m,n,mnl,f,s,er
c       m: number of grid-points
c       n: lenth of time series
c       mnl=min(m,n)
c       f(m,n): the raw spatial-temporal seires
c       s(mnl,mnl): eigenvectors
c       er(mnl,1): lamda (eigenvalues), its equence is from big to small.
c       er(mnl,2): accumulated eigenvalues from big to small
c       er(mnl,3): explained variances (lamda/total explain) from big to small
c       er(mnl,4): accumulated explaned variances from big to small
c     output: egvt,ecof
c       egvt(m,mnl): normalized eigenvectors
c       ecof(mnl,n): time coefficients of egvt 
      subroutine tcoeff(m,n,mnl,f,s,er,egvt,ecof)
      dimension f(m,n),s(mnl,mnl),er(mnl,4),egvt(m,mnl),ecof(mnl,n)
      dimension v(mnl)  !work array
      do j=1,mnl
        do i=1,m
          egvt(i,j)=0.
        enddo
        do i=1,n
          ecof(j,i)=0.
        enddo
      enddo
c-----Normalizing the input eignvectors s
      do 10 j=1,mnl
        c=0.
        do i=1,mnl
          c=c+s(i,j)*s(i,j)
        enddo
        c=sqrt(c)
        do i=1,mnl
          s(i,j)=s(i,j)/c
        enddo
  10  continue
c-----
      if(m.le.n) then
        do js=1,mnl
        do i=1,m
          egvt(i,js)=s(i,js)
        enddo
        enddo
        do 30 j=1,n
          do i=1,m
            v(i)=f(i,j)
          enddo
          do is=1,mnl
          do i=1,m
            ecof(is,j)=ecof(is,j)+v(i)*s(i,is)
          enddo
          enddo
  30    continue
      else
        do 40 i=1,m
          do j=1,n
            v(j)=f(i,j)
          enddo
          do js=1,mnl
          do j=1,n
            egvt(i,js)=egvt(i,js)+v(j)*s(j,js)
          enddo
          enddo
  40    continue
        do 50 js=1,mnl
          do j=1,n
            ecof(js,j)=s(j,js)*sqrt(abs(er(js,1)))
          enddo
          do i=1,m
            egvt(i,js)=egvt(i,js)/sqrt(abs(er(js,1)))
          enddo
  50    continue
      endif
      return
      end
c-----*----------------------------------------------------6---------7--
c     Computing the mean ax, standard deviation sx
c       and variance vx of a series x(i) (i=1,...,n).
c     input: n and x(n)
c       n: number of raw series
c       x(n): raw series
c     output: ax, sx and vx
c       ax: the mean value of x(n)
c       sx: the standard deviation of x(n)
c       vx: the variance of x(n)
c     By Dr. LI Jianping, May 6, 1998.
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