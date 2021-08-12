c11c W flux is defined by Takaya and Nakamura (1997, GRL)
        PARAMETER (mg=72,jg2=37,nl1=7,nl2=5,nm=18, L=72,M=15)
        real um(L,M),vm(L,M)
        real x(L,M),zm(72,37),z500(72,37)
        real psi_x(L,M),psi_y(L,M),ex(L,M),ey(L,M)
        real psi_xx(L,M),psi_xy(L,M),psi_yy(L,M)
        REAL dx(m),phi(m),psi(L,M),
     & cor(m),rcos(m),rsin(m),r2sin(m)
        DATA PHI1,PHI2,DELAM/20.,90.,5.0/


!! fort.41 is the climatology of a field, say Z500
!! fort.42 is an anomaly z500 pattern (regression pattern, composite pattern)

      open(25,file='fort.41',form='unformatted')
      open(26,file='fort.42',form='unformatted')
      open(23,file='Wvector.dat',form='unformatted')


CC constants
      G=9.81
      PI=4.0*ATAN(1.0)
      EARTH=6370949.0
      CRAD=PI/180.0
      DLAM=DELAM
      DPHI=(-PHI2+PHI1)/(M-1)
      RDLAM=CRAD*DLAM
      RDPHI=CRAD*DPHI
      DO 11 J=1,M
      PHI(J)=PHI2+(J-1)*DPHI
      PHI(J)=PHI(J)*CRAD
      RCOS(J)=COS(PHI(J))
      RSIN(J)=SIN(PHI(J))
        r2SIN(J)=SIN(2.*PHI(J))
      DX(J)=RDLAM*EARTH*RCOS(J)
      COR(J)=2.0*7.292E-5*RSIN(J)
   11 CONTINUE
	DY=RDPHI*EARTH

CCCCCCCC
         read(25)zm


c     Calculate geostrophic um & vm from zm
c
      do 72 j = 1,m
        f0=cor(j)
        if(j.le.15.and.j.ge.12)f0=cor(13)
      jp1 = j + 1
      jm1 = j - 1
      if (j.eq.1) jm1 = 1
      if (j.eq.m) jp1 = m
      dyy = 2.*dy
      if (j.eq.1.or.j.eq.m) dyy = dy
      do 72 i = 1,l
      ip1 = i + 1
      if (i.eq.l) ip1 = 1
      imn1 = i - 1
      if (i.eq.1) imn1 = l
      um(i,j) =-g/f0*(zm(i,jp1)-zm(i,jm1))/dyy
      vm(i,j) = g/f0*(zm(ip1,j)-zm(imn1,j))/(2.*dx(j))
   72 continue

CCCCCCCCcc--------------------------

         read(26)z500

cc get eddy streamfunction from geopotential height by linear balance equation
cc      Psi = g/f Z     f=f0 as is at 45N

          do j=1,M
           f0=cor(j)
           if(j.le.15.and.j.ge.12)f0=cor(13)
          do i=1,L
           psi(i,j)=z500(i,j)*G/f0
          end do
          end do

ccccc Psi_x, Psi_y
c
      do j = 1,m
      jp1 = j + 1
      jm1 = j - 1
      if (j.eq.1) jm1 = 1
      if (j.eq.m) jp1 = m
      dyy = 2.*dy
      if (j.eq.1.or.j.eq.m) dyy = dy
      do i = 1,l
      ip1 = i + 1
      if (i.eq.l) ip1 = 1
      imn1 = i - 1
      if (i.eq.1) imn1 = l
      psi_x(i,j) = (psi(ip1,j)-psi(imn1,j))/(2.*dx(j))
      psi_y(i,j) = (psi(i,jp1)-psi(i,jm1))/dyy
        end do
        end do


cccc Psi_xx, Psi_xy, Psi_yy
        do j = 1,m
      jp1 = j + 1
      jm1 = j - 1
      if (j.eq.1) jm1 = 1
      if (j.eq.m) jp1 = m
      dyy = 2.*dy
      if (j.eq.1.or.j.eq.m) dyy = dy
      do i = 1,l
      ip1 = i + 1
      if (i.eq.l) ip1 = 1
      imn1 = i - 1
      if (i.eq.1) imn1 = l
      psi_xx(i,j) = (psi_x(ip1,j)-psi_x(imn1,j))/(2.*dx(j))
      psi_xy(i,j) = (psi_x(i,jp1)-psi_x(i,jm1))/dyy
      psi_yy(i,j) = (psi_y(i,jp1)-psi_y(i,jm1))/dyy
        end do
        end do


cc wave-activity flux Ex, Ey
        do j = 1,m
        do i = 1,l
         absu=sqrt(um(i,j)**2+vm(i,j)**2)
         bb=1./(2.*absu)
         ex(i,j)=um(i,j)*(psi_x(i,j)**2-psi(i,j)*psi_xx(i,j))+
     *           vm(i,j)*(psi_x(i,j)*psi_y(i,j)-psi(i,j)*psi_xy(i,j))
         ey(i,j)=um(i,j)*(psi_x(i,j)*psi_y(i,j)-psi(i,j)*psi_xy(i,j))+
     *           vm(i,j)*(psi_y(i,j)**2-psi(i,j)*psi_yy(i,j))
         ex(i,j)=bb*ex(i,j)
         ey(i,j)=bb*ey(i,j)
         r=sqrt(ex(i,j)**2+ey(i,j)**2)
          if(r.lt.0.1)then
          ex(i,j)=0.
          ey(i,j)=0.
          end if
        end do
        end do



        write(23)((ex(i,j),i=1,L),j=1,M)
        write(23)((ey(i,j),i=1,L),j=1,M)





        stop
        end









