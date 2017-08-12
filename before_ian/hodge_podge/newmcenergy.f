        subroutine ljenergy(xcoor,ycoor,zcoor,xi,yi,zi,natoms,ti,boxl)

c		tgr,energy,virial=mcenergy.ljenergy(xcoor,ycoor,zcoor,natoms,boxlength,deltag,nbins)
c  v12old,v6old,w12old,w6old=newmcenergy.ljenergy(xcoor,ycoor,zcoor,rxiold,ryiold,rziold,natoms,i,boxl)
c       f2py -c -m newmcenergy newmcenergy.f

        integer natoms,ti,j
        double precision xcoor(natoms)
        double precision ycoor(natoms)
        double precision zcoor(natoms)
        double precision boxl
        double precision rcutsq,xi,yi,zi
        double precision rxij,ryij,rzij,rijsq
        double precision sr2,sr6,vij12,vij6,v12,v6,w12,w6

cf2py	intent(in) :: xcoor,ycoor,zcoor
cf2py	intent(in) :: natoms,boxl,ti,xi,yi,zi
cf2py	intent(out) :: v12,v6,w12,w6

cf2py   depend(natoms) :: xcoor,ycoor,zcoor

cf2py   intent(hide) :: j,rcutsq
cf2py   intent(hide) :: rxij,ryij,rzij,rijsq
cf2py   intent(hide) :: sr2,sr6,vij12,vij6

        v12=0.0
        v6=0.0
        w12=0.0
        w6=0.0 
        rcutsq=(boxl/2.0)**2.0      

        do 200 j=1,natoms
                if (j .ne. ti) then
                        rxij=xi-xcoor(j)
                        ryij=yi-ycoor(j)
                        rzij=zi-zcoor(j)

                        rxij=rxij-boxl*(ANINT(rxij/boxl))
                        ryij=ryij-boxl*(ANINT(ryij/boxl))
                        rzij=rzij-boxl*(ANINT(rzij/boxl))

                        rijsq=rxij*rxij+ryij*ryij+rzij*rzij
                        
                        if (rijsq .lt. rcutsq) then
                                sr2=1.0/rijsq
                                sr6=sr2*sr2*sr2
                                vij12=sr6*sr6
                                vij6=-sr6
                                v12=v12+vij12
                                v6=v6+vij6
                                w12=w12+vij12
                                w6=w6+vij6*0.5 
                        endif
                endif        
  200   continue
        v12=4.0*v12
        v6=4.0*v6
        w12=48.0*w12/3.0
        w6=48.0*w6/3.0

C        write(6,*) 'v12 = ',v12 
C        write(6,*) 'v6 = ',v6
C        write(6,*) 'w12 = ',w12 
C        write(6,*) 'w6 = ',w6
 
        return 
        end
