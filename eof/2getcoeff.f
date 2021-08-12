CC     求EOF的时间序列
      Parameter (nm=57, ng=25*29,nx=25,ny=29)
      Real t(nx,ny),tm(nx,ny),
     &      var(ng,nm), cct(ng,ng), ttp(nx,ny,nm),
     &      eign(ng,4),c(4,nm),y1(nm)

      open(15,file='
     &  F:\Data\Processed\ecmwf\era_interim.sat.partial.1958-2014DJF',
     &  form='binary',status='old')

      open(23,file='F:\Programming\fortran\eof3\eof3\eof3\
     &era_interim.sat.partial.eof',form='unformatted')
      open(13,file='F:\Data\Processed\ecmwf\
     &era_interim.sat.1958-2014DJF_PC1.txt',form='formatted')
      open(14,file='F:\Data\Processed\ecmwf\
     &era_interim.sat.1958-2014DJF_PC1',form='binary')

	do 16 i=1,nx
	do 16 j=1,ny
16	tm(i,j)=0.

	do 61 k=1,nm
         read(15)t
	 do 12 j=1,ny
	 do 12 i=1,nx
	 ttp(i,j,k)=t(i,j)
12	 tm(i,j)=tm(i,j)+t(i,j)/float(nm)
61	continue

C--------------------------------
	do 62 k=1,nm
	ii=0
        do j=1,ny
        do i=1,nx
	 ii=ii+1
         var(ii,k)=ttp(i,j,k)-tm(i,j)
	enddo
	enddo
62	continue

	read(23)eign



        do 44 i=1,4
        do 44 n=1,nm
        c(i,n)=0.0
        do 45 m=1,ng
45      c(i,n)=c(i,n)+eign(m,i)*var(m,n)
44      continue

        do kk=1,nm
          y1(kk)=c(1,kk)
        end do

          call getstd(y1,nm,std)

        do kk=1,nm
          y1(kk)=y1(kk)/std
	  write(13,*)y1(kk)*(-1.0)
        write(14)y1(kk)*(-1.0)
	  print*,1957+kk,y1(kk)*(-1.0)
        end do

        stop	
      	end



	include "F:\Programming\fortran\subroutine\getstd.f"

