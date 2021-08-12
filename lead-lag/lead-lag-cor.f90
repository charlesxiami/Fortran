program lead_lag_correlation

implicit none
integer,parameter::nt1=36,nt3=66*12-2
real nino(nt3),coef1(nt1),r1(49),r2(49),r3(49)
integer i

open(10,file='nino3.4_1949DJF-2015SON.txt')
open(11,file='hf1980-2015jfmPREC1.txt')
open(21,file='lead-lag.dat',form='unformatted')
open(22,file='lead-lag-P1.dat',form='unformatted')
open(23,file='lead-lag-P2.dat',form='unformatted')

   
read(10,*)nino
do i=1,nt1
   read(11,*)coef1(i)
end do

do i=1,49               !lead24-lag24
   call co(coef1(1:19),nino((i+28*12+1):(i+28*12+1+(19-1)*12):12),19,r2(i))
end do
do i=1,33
   call co(coef1,nino((i+28*12+1):(i+28*12+1+(nt1-1)*12):12),nt1,r1(i))
   call co(coef1(20:nt1),nino((i+47*12+1):(i+47*12+1+(17-1)*12):12),17,r3(i))
end do
do i=34,49
   r1(i)=1.0/0
   r3(i)=1.0/0
end do

write(21)r1
write(22)r2
write(23)r3
print*,nino

end program


!**********************************************************
subroutine co(x,y,t,r)
implicit none
integer t,n
real x(t),y(t),r,sxy,sx,sy,ax,ay

ax=0.0
ay=0.0
do n=1,t
   ax=ax+x(n)/real(t)
   ay=ay+y(n)/real(t)
end do

sxy=0.0
sx=0.0
sy=0.0
do n=1,t
   sxy=sxy+x(n)*y(n)
   sx=sx+(x(n)-ax)**2
   sy=sy+(y(n)-ay)**2
end do
r=(sxy-t*ax*ay)/sqrt(sx*sy)
return

end subroutine


