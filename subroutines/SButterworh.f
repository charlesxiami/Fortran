      subroutine Bfilter2(n,x,y,a,b1,b2)
!C     Second Order Butterworth Band-Pass Filter
!C         y(k)=S(j=0,N)a(j)*x(k-j)-S(j:1,L)b(j)*y(k-j)
!C       where L=2 is called the order of the filter, and generally, N=L.
!C
!C **  Note that: Call the subroutine responseBF2 first before calling this subroutine.
!C
!C     Input parameters: n, x(n), a, b1, and b2
!C        n: the number od data
!C        x: the original series
!C        a, b1 and b2 from the subroutine responseBF2
!C     Output variable: y(n)
!C        y: the final filtered series of x
!C     Work array: y1(n)
!C        y1: the initial filtered series of x, work array 
!c     By Dr. LI Jianping, March 8, 2000.
!c-----*----------------------------------------------------6---------7--
      dimension x(n),y(n)
      dimension y1(n)                                !Work array
      y1(1)=0.
      y1(2)=0.
      do 10 i=3,n
        y1(i)=a*(x(i)-x(i-2))-(b1*y1(i-1)+b2*y1(i-2))
  10  continue
      y(n)=y1(n)
      y(n-1)=y1(n-1)
      do 20 i=n-2,1,-1
        y(i)=a*(y1(i)-y1(i+2))-(b1*y(i+1)+b2*y(i+2))
  20  continue
  30  continue
      return
      end
!c-----*----------------------------------------------------6---------7--
      subroutine responseBF2(n,dt,w1,w2,a,b1,b2,h)
!C     Frequency Response Function of Second-order Butterworth Band-pass Filter
!C     
!C **  Note that: Call the subroutine responseBF2 first before call the subroutine
!C                Bfilter2() for the Second Order Butterworth Band-Pass Filter
!C
!C     Input parameters: n, dt, w1, and w2
!C       dt: the sampling interval (e.g., dt=1. if you use daily data)
!C       w1: lower (circular) cutoff frequency on the left of w0
!C           w1=2*pi/t1 where t1 is corresponding period to w1
!C       w2: upper (circular) frequency on the right of w0
!C           w2=2*pi/t2 where t2 is corresponding period to w2 
!C       w0: central (circular) frequency
!C           w0=sqrt(w1*w2), =2*pi/t0 where t0 is corresponding period to w0 
!C       w1<w0<w2, that is t1>t0>t2
!C     Output variables: a, b1, b2, and h
!C        h: =|w(z)|**2 the amplitude of the frequency response function w(z)
!C           where z=exp(-i*omg*datt). Note that: For F90, z=exp((0,-omg)) should 
!C           be modified to F90 format, z=exp(cmplx(0,-omg)). 
!C           Here we compile it in F77 format.
!C     Work array:
!C        w: the frequency response function w(z), it is a complex array. 
!C     By Dr. LI Jianping, March 8, 2000. 
      dimension h(n)
      complex w,z                       !Work array
      pi=3.1415926
      w0=sqrt(w1*w2)
      dt=1.
      a1=sin(w1*dt)/(1.+cos(w1*dt))
      a2=sin(w2*dt)/(1.+cos(w2*dt))
      dw=2.*abs(a1-a2)
      omg2=4.*a1*a2
      c=4.+2.*dw+omg2
      a=2.*dw/c
      b1=2.*(omg2-4.)/c
      b2=(4.-2.*dw+omg2)/c
      do 10 i=1,n
        omg=2.*pi/float(i)
        omg=omg*dt
!        z=exp((0,-omg)) !F77
        z=exp(cmplx(0,-omg)) !F90
        w=a*(1-z**2)/(1.+b1*z+b2*z**2)
        h(i)=abs(w)**2
  10  continue
      return
      end