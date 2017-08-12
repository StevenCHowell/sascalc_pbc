        subroutine pairsans(xcoor,ycoor,zcoor,qx,qy,
     +qz,batom,natoms,numqpoints,nqvectors,sq)

        double precision xcoor(natoms)
        double precision ycoor(natoms)
        double precision zcoor(natoms)
        double precision qx(numqpoints*nqvectors)
        double precision qy(numqpoints*nqvectors)
        double precision qz(numqpoints*nqvectors)
        double precision batom(natoms)
        double precision sq(numqpoints)
        integer natoms,numqpoints,nqvectors

        double precision sum,tcosu,tcosv,tsinu,tsinv
        double precision tqx,tqy,tqz,b1,b2,val
        integer q,qvec,i,j,k,qs        

cf2py	intent(in) :: xcoor,ycoor,zcoor
cf2py	intent(in) :: qx,qy,qz,batom
cf2py	intent(out) :: sq

cf2py   depend(natoms) :: xcoor,ycoor,zcoor,batom
cf2py   depend(numqpoints,nqvectors) :: qx,qy,qz
cf2py   depend(numqpoints) :: sq

cf2py   intent(hide) :: sum,tcosu,tcosv,tsinu,tsinv
cf2py   intent(hide) :: q,qvec,i,j,k,qs
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
               u=tqx*xcoor(i)+tqy*ycoor(i)+tqz*zcoor(i)
               tcosu=cos(u)
               tsinu=sin(u)
               do 100,j=1,natoms
                 b2=batom(j)
                 v=tqx*xcoor(j)+tqy*ycoor(j)+tqz*zcoor(j)
                 tcosv=cos(v)
                 tsinv=sin(v)
                 val=b1*b2*(tcosu*tcosv+tsinu*tsinv)
                 sum=sum+val

  100          continue
  200        continue
             qs=qs+1 
  300      continue
c           sq(q)=sum/float(nqvectors*natoms)
           sq(q)=sum/float(nqvectors)
  400   continue
        return 
        end
