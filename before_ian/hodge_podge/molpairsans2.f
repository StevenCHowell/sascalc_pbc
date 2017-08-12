        subroutine pairsans(xcoor,ycoor,zcoor,qx,qy,
     +qz,batom,natoms,numqpoints,nqvectors,boxl,
     +sq)

        double precision xcoor(natoms)
        double precision ycoor(natoms)
        double precision zcoor(natoms)
        double precision qx(numqpoints*nqvectors)
        double precision qy(numqpoints*nqvectors)
        double precision qz(numqpoints*nqvectors)
        double precision batom(natoms)
        double precision boxl
        double precision sq(numqpoints)
        complex arg
        integer natoms,numqpoints,nqvectors

        double precision sum,tqx,tqy,tqz,b1,b2,val
        integer q,qvec,i,j,k,qs        

cf2py	intent(in) :: xcoor,ycoor,zcoor
cf2py	intent(in) :: qx,qy,qz,batom,boxl
cf2py	intent(out) :: sq

cf2py   depend(natoms) :: xcoor,ycoor,zcoor,batom
cf2py   depend(numqpoints,nqvectors) :: qx,qy,qz
cf2py   depend(numqpoints) :: sq

cf2py   intent(hide) :: sum,q,qvec,i,j,k,qs,arg
cf2py   intent(hide) :: tqx,tqy,tqz,b1,v2,val

        do 10,q=1,numqpoints
           sq(q)=0.0
  10    continue

C        do 20,i=1,nqvectors*numqpoints
C           write(*,*) qx(i),qy(i),qz(i)
C  20    continue
        k=1
        qs=1
        do 400,q=1,numqpoints
           write(*,*) 'q = ',q
           sum=0.0
           val=0.0
           do 300,qvec=1,nqvectors
             tqx=qx(qs)
             tqy=qy(qs)
             tqz=qz(qs)
             do 200,i=1,natoms
               b1=batom(i)
               xi=xcoor(i)
               yi=ycoor(i)
               zi=zcoor(i)
               do 100,j=1,natoms
                 b2=batom(j)
                 xj=xcoor(j)
                 yj=ycoor(j)
                 zj=zcoor(j)
                 rxij=xi-xj
                 ryij=yi-yj
                 rzij=zi-zj
                 rxij=rxij-boxl*(ANINT(rxij/boxl))
                 ryij=ryij-boxl*(ANINT(ryij/boxl))
                 rzij=rzij-boxl*(ANINT(rzij/boxl))
                 arg=EXP((0.0,-1.0)*(tqx*rxij+tqy*ryij+tqz*rzij))
                 if (j .ne. i) then
                   val=b1*b2*REAL(arg)
                   sum=sum+val
                 endif

  100          continue
  200        continue
             qs=qs+1 
  300      continue
c           sq(q)=sum/float(nqvectors*natoms)
           sq(q)=sum/float(nqvectors)
  400   continue
        return 
        end
