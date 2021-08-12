subroutine runningmean(x,y,n,a)
integer n,a
real::y(n),x(n)
real s
integer i,j,k
!*********************************************
!x: array of index
!y: running mean of x
!n: samlpe size
!a: sliding window

do i=1,(a-1)/2
   do j=i,n-i+1,n-i*2+1
      y(j)=0.0
      do k=j-i+1,j+i-1
         y(j)=y(j)+x(k)/real(i*2-1)
	  end do
   end do
end do

i=(a+1)/2
do j=i,n-i+1
   y(j)=0.0
   do k=j-i+1,j+i-1
      y(j)=y(j)+x(k)/real(a)
   end do
end do


return
end subroutine runningmean