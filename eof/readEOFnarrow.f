        parameter(lon=25,lat=25,ng=lon*lat, nMax=4)
        real x(lon,lat),evt(ng,nMax)


      open(17,file='EOF_phi.dat',form='binary')
      read(17) evt
      open(30,file='EOF_phi_pattern.dat',form='binary')
         do 99 k=1,4
         ii=0
         do 14 j=1,lat
         do 14 i=1,lon
          ii=ii+1
          x(i,j)=evt(ii,k)
14      continue
         write(30)((x(i,j),i=1,lon),j=1,lat)
         
99	 continue

	end
