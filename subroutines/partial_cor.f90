!get partial regression
subroutine partial_cor(x1,x2,y,n,pa,pr)
integer n,i
real x1(n),x2(n),y(n),pr,xres(n)
real(8) xmean,ymean,sxy,sx,sy,r,pa
!*********************************************
!x1: array of the first index
!x2: array of the second index
!remove the second index from the first index
!y: array of index or pattern
!n: samlpe size
!pa: coefficient of  partial regression
!pa: partial correlation coefficient

call correlation(x1,x2,n,r)

do i=1,n
   xres(i)=x1(i)-r*x2(i)
end do

call correlation(xres,y,n,pr)
call regression(xres,y,n,pa)

return
end subroutine partial_cor



include "correlation.f90"
include "regression.f90"