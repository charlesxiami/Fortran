subroutine ttest(x1,x2,n,m,t)
implicit none
integer n,m
real x1(n),x2(m)
real mean1,mean2,sx1,sx2,s,t,l
integer i
mean1=0.0
mean2=0.0
do i=1,n
    mean1=mean1+x1(i)/float(n)
end do
do i=1,m
    mean2=mean2+x2(i)/float(m)
end do
!print*,mean1,mean2
sx1=0.0
sx2=0.0
do i=1,n
    sx1=sx1+(x1(i)-mean1)*(x1(i)-mean1)
end do
do i=1,m
    sx2=sx2+(x2(i)-mean2)*(x2(i)-mean2)
end do

s=sqrt((sx1+sx2)/(real(n+m-2)))
t=(mean1-mean2)/(s*sqrt(1.0/n+1.0/m))

return
end subroutine ttest