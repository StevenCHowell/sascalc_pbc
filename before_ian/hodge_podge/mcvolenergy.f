        subroutine ljvenergy(xcoor,ycoor,zcoor,
     +natoms,boxl,energy)

c		energy=mcvolenergy.ljvenergy(xcoor,ycoor,zcoor,natoms,boxlength)
c       f2py -c -m mcvolenergy mcvolenergy.f

        double precision xcoor(natoms)
        double precision ycoor(natoms)
        double precision zcoor(natoms)
        double precision boxl
        integer natoms

        double precision val,xi,yi,zi,xj,yj,zj
        double precision rxij,ryij,rzij,rij
        double precision irij2,irij6,irij12
        double precision corr,rij2,density
        double precision irij3,irij9
        double precision energy

        integer i,j

cf2py	intent(in) :: xcoor,ycoor,zcoor
cf2py	intent(in) :: natoms,boxl
cf2py	intent(out) :: energy

cf2py   depend(natoms) :: xcoor,ycoor,zcoor

cf2py   intent(hide) :: val,i,j,k,xi,yi,zi
cf2py   intent(hide) :: rxij,ryij,rzij,rij
cf2py   intent(hide) :: irij2,irij6,irij12
cf2py   intent(hide) :: rij2,irij3,irij9
cf2py   intent(hide) :: corr,density

        energy=0.0
        density=real(natoms)/(boxl**3.0)
        corr=8.0*3.1415926*density/3.0 
        do 200,i=1,natoms-1
                xi=xcoor(i)
                yi=ycoor(i)
                zi=zcoor(i)
                do 100,j=i+1,natoms
                        xj=xcoor(j)
                        yj=ycoor(j)
                        zj=zcoor(j)
                        rxij=xi-xj
                        ryij=yi-yj
                        rzij=zi-zj
                        rxij=rxij-boxl*(ANINT(rxij/boxl))
                        ryij=ryij-boxl*(ANINT(ryij/boxl))
                        rzij=rzij-boxl*(ANINT(rzij/boxl))

                        rij=SQRT(rxij*rxij+ryij*ryij+rzij*rzij)
                        rij2=rij*rij
                        if(rij .lt. boxl/2.0) then
                                irij2=1.0/rij2
                                irij6=irij2**3.0
                                irij12=irij6**2.0
                                val=4.0*(irij12-irij6) 
                                
                        else
                                irij=1.0/rij
                                irij3=irij**3.0
                                irij9=irij3**3.0
                                val=corr*((irij9/3.0)-irij3) 
                        endif
                        energy=energy+val
  100           continue
  200   continue
        
        return 
        end
