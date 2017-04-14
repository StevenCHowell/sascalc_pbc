c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine calc_gr(coor,boxl,gr,dr,natoms,nbins)
        double precision boxl,dr
        double precision coor(natoms,3)
        integer natoms,nbins,i,j,ri
        double precision gr(nbins)
        double precision dxij,dyij,dzij,rijsq,rcutsq

cf2py	intent(in):: coor,boxl,dr,natoms,nbins
cf2py	intent(in):: gr
cf2py	intent(hide):: i,j,ri
cf2py	intent(hide):: dxij,dyij,dzij,rijsq,rcutsq

        rcutsq=(boxl/2.0)*(boxl/2.0)

        do 200,i=1,natoms-1
                do 300,j=i+1,natoms
                        dxij=coor(i,1)-coor(j,1)
                        dyij=coor(i,2)-coor(j,2)
                        dzij=coor(i,3)-coor(j,3)

                        dxij=dxij-boxl*(ANINT(dxij/boxl))
                        dyij=dyij-boxl*(ANINT(dyij/boxl))
                        dzij=dzij-boxl*(ANINT(dzij/boxl))

                        rijsq=dxij*dxij+dyij*dyij+dzij*dzij

                        if (rijsq .lt. rcutsq) then
                           ri = int(dsqrt(rijsq)/dr) + 1
                           gr(ri) = gr(ri) + 2.0
                        endif

  300   continue
  200   continue

        end

c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012
c
c     to setup and incorporate into python:
c
c     python setup_gr.py build
c     cp build/lib*/distance.so ./
c
c     to call this from python:
c
c     import sys ; sys.path.append('./')
c     from distance import distance
c
c     dist = distance(coor, dist)
c
c         1         2         3         4         5         6         7
c123456789012345678901234567890123456789012345678901234567890123456789012
