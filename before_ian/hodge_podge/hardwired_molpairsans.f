        subroutine pairsans(xcoor,ycoor,zcoor,natoms)
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

        double precision sum,cosu,cosv,sinu,sinv
        double precision tqx,tqy,tqz,b1,b2,val
        integer q,qvec,i,j,k,qs        

cf2py	intent(in) :: xcoor,ycoor,zcoor
cf2py	intent(in) :: qx,qy,qz,batom
cf2py   intent(in) :: natoms,numqpoints,nqvectors
cf2py	intent(out) :: sq
cf2py   intent(hide) :: sum,cosu,cosv,sinu,sinv
cf2py   intent(hide) :: q,qvec,i,j,k,qs
cf2py   intent(hide) :: tqx,tqy,tqz,b1,v2,val

        do 10,q=1,numqpoints
           sq(q)=0.0
  10    continue

        k=0
        qs=0
        do 400,q=1,numqpoints
           sum=0.0
           do 300,qvec=1,nqvectors
             tqx=qx(qs)
             tqy=qy(qs)
             tqz=qz(qs)
             do 200,i=1,natoms
               b1=batom(i)
               u=tqx*xcoor(i)+tqy*ycoor(i)+tqz*zcoor(i)
               cosu=cos(u)
               sinu=sin(u)
               do 100,j=1,natoms
                 b2=batom(j)
                 v=tqx*xcoor(j)+tqy*ycoor(j)+tqz*zcoor(j)
                 cosv=cos(v)
                 sinv=sin(v)
                 val=b1*b2*(cosu*cosv+sinu*sinv)
                 sum=sum+val

  100          continue
  200        continue
             qs=qs+1 
  300      continue
           sq(k)=sum/float(nqvectors*natoms)
  400   continue
        return 
        end
