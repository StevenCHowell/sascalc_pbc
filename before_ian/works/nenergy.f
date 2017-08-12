
        subroutine lje(coor,x1,y1,z1,i,boxl,natoms,deltag,nbins,
     +tgr,lv12,lv6,lw12,lw6)
        double precision coor(natoms,3)
        double precision boxl
        double precision x1,y1,z1,lv12,lv6,lw12,lw6
        integer natoms,i,j,nbins,ig
        double precision tgr(nbins)
        double precision x2,y2,z2
        double precision sr2,sr6,vij12,vij6,rcutsq
        double precision rxij,ryij,rzij,rijsq,deltag

cf2py	intent(in)  coor,x1,y1,z1,num,boxl,deltag,nbins
cf2py	intent(out) tgr,lv12,lv6,lw12,lw6
cf2py	intent(hide) natoms,j,ig
cf2py   depend(nbins) tgr
cf2py	intent(hide) x2,y2,z2,rxij,ryij,rzij,rijsq,rcutsq

        lv12=0.0
        lv6=0.0
        lw12=0.0
        lw6=0.0
        rcutsq=(boxl/2.0)**2.0

        do 200,j=1,natoms
                if(j .ne. i) then 
                        x2=coor(j,1)
                        y2=coor(j,2)
                        z2=coor(j,3)
         
                        rxij=x2-x1
                        ryij=y2-y1
                        rzij=z2-z1
         
                        rxij=rxij-boxl*(ANINT(rxij/boxl))
                        ryij=ryij-boxl*(ANINT(ryij/boxl))
                        rzij=rzij-boxl*(ANINT(rzij/boxl))
                         
                        rijsq=rxij*rxij+ryij*ryij+rzij*rzij

                        if (rijsq .lt. rcutsq) then
                                sr2=1.0/rijsq
                                sr6=sr2*sr2*sr2
                                vij12=sr6*sr6
                                vij6=-sr6
                                lv12=lv12+vij12
                                lv6=lv6+vij6
                                lw12=lw12+vij12
                                lw6=lw6+vij6*0.5
                                ig = int(sqrt(rijsq)/deltag)
                                tgr(ig) = tgr(ig) + 2.0
                        endif
                endif
  200   continue

        lv12=4.0*lv12
        lv6=4.0*lv6
        lw12=48.0*lw12/3.0
        lw6=48.0*lw6/3.0

        end
  
