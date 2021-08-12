!        EOF：经验正交函数分解!                “实现时空场的分离”具体可以参考黄嘉佑老师的《气象统计分析与预报方法》
!        这是来自“lijianping”老师的EOF分解程序，在此对李建平老师表示感谢O(∩_∩)O~
!        相关参数说明
!        m：格点数
!        n：时间长度
!        mnl：模态数                        
!        ks：[1]-1表示对原场EOF；
!                [2]0表示多距平场EOF；
!                [3]1表示对归一化场EOF。
!        "test.txt"        ：需要分解的时空场
!        "er.txt"：存放解释方差        [1]：er(mnl,1)特征值——        [2]：er(mnl,2)累计特征值——        [3]：er(mnl,3)解释方差——        [4]：er(mnl,4)累计解释方差!【20120922补充：解释方差就是对特征值的归一化】
!        "ecof.txt"：存放时间系数------每一列是一个模态
!        "egvt.txt"：存放各个模态的空间向量场------每一列是一个模态
program main 

integer, parameter:: m=160,n=61,mnl=61,ks=1

real x(m,n),egvt(m,mnl),ecof(mnl,n),er(mnl,4)
 open(1,file='e:/eof/bzh.txt',form='formatted')
 read(1,100)((x(i,j),i=1,m),j=1,n)
 close(1)
call eof(m,n,mnl,x,ks,er,egvt,ecof)
        open(2,file='e:/eof/er.txt')
        do i=1,mnl
        write(2,"(4f30.8)") (er(i,j),j=1,4)
        enddo
        close(2)
        open(3,file='e:/eof/ecof.txt')
        do j=1,n
        write(3,"(<mnl>f20.6)") (ecof(i,j),i=1,mnl)
        enddo
        close(3)
        open(4,file='e:/eof/egvt.txt')
        do i=1,m
        write(4,"(<mnl>f20.6)") (egvt(i,j),j=1,mnl)
        enddo
        close(4)
print*,"-------------------------------"
print*,"这是来自李建平老师的EOF分解程序"
print*,"http://bbs.06climate.com整理制作"
print*,"-------------------------------"
100 format(5f15.7)
end

!***************************************************************************
        subroutine eof(m,n,mnl,f,ks,er,egvt,ecof)
                implicit none
                integer*4                                ::        m,n,mnl,ks
                real*4                                        ::        f(m,n),er(mnl,4),egvt(m,mnl),ecof(mnl,n)
                real*4,allocatable                ::        cov(:,:),s(:,:),d(:)
                call transf(m,n,f,ks)
                allocate(cov(mnl,mnl))
                call crossproduct(m,n,mnl,f,cov)
                allocate(s(mnl,mnl))
                allocate(d(mnl))
                call jacobi(mnl,cov,s,d,0.00001)
                call arrang(mnl,d,s,er)
                call tcoeff(m,n,mnl,f,s,er,egvt,ecof)
                return
        end
        subroutine transf(m,n,f,ks)
                implicit none
                integer*4::        m,n,ks
                real*4 ::        f(m,n)
                real*4,allocatable  ::        fw(:),wn(:)
                integer*4   ::        i0,i,j
                real*4   ::        af,sf,vf
                allocate(fw(n))
                allocate(wn(m))
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
                        write(*,*)'*************************  fault  ***********************************'
                        write(*,*)'The program cannot go on because the original field has invalid data.'
                        write(*,*)'There are totally ',i0,'gridpionts with invalid data.'
                        write(*,*)'The array wn(m) stores the positions of those invalid grid-points.'
                        write(*,*)'You must pick off those invalid data from the orignal field---$ '  
                        write(*,*)'$---and then reinput a new field to calculate its eofs.' 
                        write(*,*)'************************  fault  ************************************'
                endif            
                if(ks.eq.-1)return
                if(ks.eq.0)then
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
                if(ks.eq.1)then
                        do i=1,m
                        do j=1,n
                                fw(j)=f(i,j)
                        enddo
                        call meanvar(n,fw,af,sf,vf)
                        if(sf==0.0)then
                        do j=1,n
                        f(i,j)=0.0
                        enddo
                        else
                        do j=1,n
                        f(i,j)=(f(i,j)-af)/sf
                        enddo
                        endif
                        enddo
                endif
                return
        end
        subroutine crossproduct(m,n,mnl,f,cov)
                implicit none
                integer*4                ::        m,n,mnl
                real*4                        ::        f(m,n),cov(mnl,mnl)
                integer*4                ::        i,j,is,js
                        if((n-m)<0)then
                        do i=1,mnl
                        do j=i,mnl
                                cov(i,j)=0.0
                                do is=1,m
                                        cov(i,j)=cov(i,j)+f(is,i)*f(is,j)
                                enddo
                                cov(j,i)=cov(i,j)
                        enddo
                        enddo
                        else
                        do i=1,mnl
                        do j=i,mnl
                                cov(i,j)=0.0
                                do js=1,n
                                        cov(i,j)=cov(i,j)+f(i,js)*f(j,js)
                                enddo
                                cov(j,i)=cov(i,j)
                        enddo
                        enddo
                        endif
                return        
        end
        subroutine jacobi(m,a,s,d,eps)
                implicit none
                integer*4                ::        m
                real*4                        ::        a(m,m),s(m,m),d(m),eps
                integer*4                ::        i,j,i1,l,iq,iq1,ip
                real*4                        ::        g,s1,s2,s3,v1,v2,v3,u,st,ct
                        s=0.0
                        do i=1,m
                                s(i,i)=1.0
                        enddo
                        g=0.
                        do i=2,m
                                i1=i-1
                                do j=1,i1
                                g=g+2.0*a(i,j)*a(i,j)
                                enddo
                        enddo
                        s1=sqrt(g)
                        s2=eps/float(m)*s1
                        s3=s1
                        l=0
50                        s3=s3/float(m)
60                        do iq=2,m
                                iq1=iq-1                                        
                                do ip=1,iq1
                                        if(abs(a(ip,iq)).lt.s3) exit
                                        l=1
                                        v1=a(ip,ip)
                                        v2=a(ip,iq)
                                        v3=a(iq,iq)
                                        u=0.5*(v1-v3)
                                        if(u.eq.0.0) g=1.
                                        if(abs(u).ge.1e-10) g=-sign(1.,u)*v2/sqrt(v2*v2+u*u)
                                        st=g/sqrt(2.*(1.+sqrt(1.-g*g)))
                                        ct=sqrt(1.-st*st)
                                        do i=1,m
                                                g=a(i,ip)*ct-a(i,iq)*st
                                                a(i,iq)=a(i,ip)*st+a(i,iq)*ct
                                                a(i,ip)=g
                                                g=s(i,ip)*ct-s(i,iq)*st
                                                s(i,iq)=s(i,ip)*st+s(i,iq)*ct
                                                s(i,ip)=g
                                        enddo
                                        do i=1,m
                                                a(ip,i)=a(i,ip)
                                                a(iq,i)=a(i,iq)
                                        enddo
                                        g=2.*v2*st*ct
                                        a(ip,ip)=v1*ct*ct+v3*st*st-g
                                        a(iq,iq)=v1*st*st+v3*ct*ct+g
                                        a(ip,iq)=(v1-v3)*st*ct+v2*(ct*ct-st*st)
                                        a(iq,ip)=a(ip,iq)
                                enddo
                        enddo
                        if((l-1)==0)then
                                l=0
                                go to 60
                        else
                                go to 150
                        endif
150                        if(s3.gt.s2) then 
                         goto 50
                        end if
                        do i=1,m
                                d(i)=a(i,i)
                        enddo
                return
        end
        subroutine arrang(mnl,d,s,er)
                implicit none
                integer*4                ::        mnl
                real*4                        ::        d(mnl),s(mnl,mnl),er(mnl,4)
                integer*4                ::        i,mnl1,k1,k2
                real*4                        ::        c,tr
                tr=0.0
                do i=1,mnl
                        tr=tr+d(i)
                        er(i,1)=d(i)
                enddo
                mnl1=mnl-1
                do k1=mnl1,1,-1
                        do k2=k1,mnl1
                        if(er(k2,1).lt.er(k2+1,1)) then
                        c=er(k2+1,1)
                        er(k2+1,1)=er(k2,1)
                        er(k2,1)=c                        
                        do i=1,mnl
                                c=s(i,k2+1)
                                s(i,k2+1)=s(i,k2)
                                s(i,k2)=c
                        enddo
                        endif
                        enddo
                enddo
                er(1,2)=er(1,1)
                do i=2,mnl
                        er(i,2)=er(i-1,2)+er(i,1)
                enddo
                do i=1,mnl
                        er(i,3)=er(i,1)/tr
                        er(i,4)=er(i,2)/tr
                enddo
                return
        end
        subroutine tcoeff(m,n,mnl,f,s,er,egvt,ecof)
                implicit none
                integer*4                                        ::        m,n,mnl
                real*4                                                ::        f(m,n),s(mnl,mnl),er(mnl,4),egvt(m,mnl),ecof(mnl,n)
                real*4,allocatable                        ::        v(:)
                integer*4                                        ::        i,j,js,is
                real*4                                                ::        c
                allocate(v(mnl))
                do j=1,mnl
                        do i=1,m
                                egvt(i,j)=0.
                        enddo
                        do i=1,n
                                ecof(j,i)=0.
                        enddo
                enddo
                do j=1,mnl
                        c=0.
                        do i=1,mnl
                                c=c+s(i,j)*s(i,j)
                        enddo
                        c=sqrt(c)
                        do i=1,mnl
                                s(i,j)=s(i,j)/c
                        enddo
                enddo
                if(m.le.n) then
                        do js=1,mnl
                                do i=1,m
                                        egvt(i,js)=s(i,js)
                                enddo
                        enddo
                        do j=1,n
                                do i=1,m
                                v(i)=f(i,j)
                                enddo
                                do is=1,mnl
                                do i=1,m
                                ecof(is,j)=ecof(is,j)+v(i)*s(i,is)
                                enddo
                                enddo
                        enddo
                else
                        do i=1,m
                                do j=1,n
                                v(j)=f(i,j)
                                enddo
                                do js=1,mnl
                                do j=1,n
                                egvt(i,js)=egvt(i,js)+v(j)*s(j,js)
                                enddo
                                enddo
                        enddo
                        do js=1,mnl
                                do j=1,n
                                ecof(js,j)=s(j,js)*sqrt(abs(er(js,1)))
                                enddo
                                do i=1,m
                                egvt(i,js)=egvt(i,js)/sqrt(abs(er(js,1)))
                                enddo
                        enddo
                endif
                return        
        end
        subroutine meanvar(n,x,ax,sx,vx)
                implicit none
                integer*4                ::        n
                real*4                        ::        x(n)
                real*4                        ::        ax,vx,sx
                integer*4                ::        i
                ax=0.
                vx=0.
                sx=0.
                do i=1,n
                        ax=ax+x(i)
                enddo                        
                ax=ax/float(n)
                do i=1,n
                        vx=vx+(x(i)-ax)**2
                enddo
                vx=vx/float(n)
                sx=sqrt(vx)
                return
        end
