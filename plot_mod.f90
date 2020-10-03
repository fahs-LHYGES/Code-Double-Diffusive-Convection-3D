MODULE plot_module
  use residual_module
contains
  subroutine get_plot()
    use maillage
    implicit none 
    real*8  :: dep ,x,y,z
    integer :: s,i,j,k

    ne=20
    dep = 1.d0/dble(ne)

    nbn=(ne+1)* (ne+1) * (ne+1)

    ALLOCATE(xyz(3,nbn))

    x = 0.d0
    s=1
    do i = 1,ne +1  
       y = 0.d0
       do j = 1,ne +1 
          z = 0.d0
          do k = 1,ne +1 
             xyz(1,s) = x         
             xyz(2,s) = y
             xyz(3,s) = z

             s= s + 1 

             z = z + dep
          end do
          y = y + dep 
       end do
       x= x + dep 
    end do

    allocate(phi_x(nbn),phi_y(nbn),vx(nbn),vy(nbn),vz(nbn),conc(nbn))
    allocate (temp(nbn),sher(nbn), Nusselt(nbn))
  end subroutine get_plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine tecplot (i)
    USE maillage
    IMPLICIT NONE
    INTEGER                             :: i


    WRITE(*,*) ' Tecplot !!'
    if (i.eq.1) then 
       WRITE(27,*)'TITLE = "analytical" '
       WRITE(27,*)'VARIABLES = "X", "Y",  "Z", "Phi_X","Phi_Y","vx","vy","vz","CK", "temp", "nul"'
       WRITE (27,*) ' ZONE T = RES1 I=', ne+1,  '  J=', ne+1,'    K=', ne+1,' DATAPACKING=POINT'
       do i=1,nbn
          write (27,*)  xyz(1,i), xyz(2,i), xyz(3,i) , phi_x(i), phi_y(i), vx(i),vy(i),vz(i),conc(i),temp(i),'0.000'
       end do
    else 
       WRITE (27,*)  'ZONE T =RES1 I=', ne+1,  '  J=', ne+1,'    K=', ne+1,' DATAPACKING=POINT'
       do i=1,nbn
          write (27,*)  xyz(1,i), xyz(2,i), xyz(3,i) , phi_x(i), phi_y(i), vx(i),vy(i),vz(i),conc(i),temp(i),'0.000'
       end do
    end if
  END Subroutine tecplot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_Phi_C_vx_vy (sss,X)
    use in_data
    use maillage
    use com_module
    implicit none
    real*8 , dimension(:)  ,intent(in)  :: X(nve+nvg)
    integer                             :: e,i,j,k,l,m,n,u,v,w,s,p,t,sss
    real*8                              :: xx,yy,zz
    real*8                              :: ri,rj,rk,rl,rm,rn,ru,rv,rw,rs,rp,rt
    real*8 , dimension(:,:)             :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn),EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt)
    real*8                              :: sher_0,sher_1,Nusselt_0,Nusselt_1



    call vect_to_Matr (X,EE,GG) 

    call Calc_AA(EE,GG,AA)
    call Calc_BB(EE,GG,BB)

    do e=1,nbn 
       xx = xyz(1,e)
       yy = xyz(2,e)
       zz = xyz(3,e)
       phi_x(e)  = 0.d0
       phi_y(e)  = 0.d0
       vx(e)   = 0.d0
       vy(e)   = 0.d0
       vz(e)   = 0.d0   !QIAN!!!!!!!!!!!!!!!!!!
       do i=0,ni 
          do j=1,nj
             do k=1,nk
                ri = dble(i)
                rj = dble(j)
                rk = dble(k)
                phi_x(e) = phi_x(e) + AA(i,j,k) * cos(ri*pi*xx) * sin(rj*pi*yy) * sin(rk*pi*zz)
                vy (e) = vy (e) + rk * AA(i,j,k) * cos(ri*pi*xx) * sin(rj*pi*yy) * cos(rk*pi*zz) 
                vz (e) = vz (e) - rj * AA(i,j,k) * cos(ri*pi*xx) * cos(rj*pi*yy) * sin(rk*pi*zz)
             end do
          end do
       end do

       do l=1,nl 
          do m=0,nm
             do n=1,nn
                rl = dble(l)
                rm = dble(m)
                rn = dble(n)
                phi_y(e) = phi_y(e) + BB(l,m,n) * sin(rl*pi*xx) * cos(rm*pi*yy) * sin(rn*pi*zz)
                vx (e) = vx (e) - rn * BB(l,m,n) * sin(rl*pi*xx) * cos(rm*pi*yy) * cos(rn*pi*zz)
                vz (e) = vz (e) + rl * BB(l,m,n) * cos(rl*pi*xx) * cos(rm*pi*yy) * sin(rn*pi*zz)
             end do
          end do
       end do

       vx(e)  = vx(e) * pi
       vy(e)  = vy(e) * pi
       vz(e)  = vz(e) * pi

       conc(e) = 0.d0
       sher(e) = 0.D0

       do u=0,nu
          do v=1,nv
             do w=0,nw
                ru = dble(u)
                rv = dble(v)
                rw = dble(w) 
                conc(e) = conc(e) + EE(u,v,w)* cos(ru*pi*xx) * sin(rv*pi*yy) * cos(rw*pi*zz)  
                if (yy.lt.1.d-4) then
                   sher(e) = sher(e) + pi*rv * EE(u,v,w) * cos(ru*pi*xx) * cos(rw*pi*zz)
                elseif (yy.gt.0.9999) then
                   sher(e) = sher(e) + pi*rv * EE(u,v,w) * cos(ru*pi*xx) * cos(rv*pi) * cos(rw*pi*zz)
                endif
             end do
          end do
       end do
       conc(e)  =  conc(e) + (1.D0-yy) 
       sher(e)  =  sher(e) - 1.d0 


       temp(e) = 0.d0
       Nusselt(e)=0.D0

       do s=1,ns
          do p=0,np
             do t=0,nt
                rs = dble(s)
                rp = dble(p)
                rt = dble(t) 
                temp(e) = temp(e) + GG(s,p,t)* sin(rs*pi*xx) * cos(rp*pi*yy) * cos(rt*pi*zz) 
                if (xx.lt.1.d-4) then
                   Nusselt(e)=Nusselt(e)+pi*rs*GG(s,p,t)*cos(rp*pi*yy)*cos(rt*pi*zz)
                elseif (yy.gt.0.9999) then
                   Nusselt(e)=nusselt(e)+pi*rs*GG(s,p,t)*cos(rs*pi)*cos(rp*pi*yy)*cos(rt*pi*zz)
                end if
             end do
          end do
       end do
       temp(e)  =  temp(e) + (1.D0-xx) 
       Nusselt(e)=Nusselt(e)-1.d0


       if (yy.lt.1.d-4) write(33,*) xx,zz,sher(e)
       if (yy.gt.0.9999) write(34,*) xx,zz,sher(e)
       if (xx.lt.1.d-4) write(35,*) yy,zz,Nusselt(e)
       if (xx.gt.0.9999) write(36,*) yy,zz,Nusselt(e)
    end do


    !Average sherwood
    sher_0=0.D0
    sher_1=0.D0
    do v=1,nv
       rv = dble(v)
       sher_0 = sher_0 + pi * rv * EE(0,v,0)
       sher_1 = sher_1 + pi * rv * EE(0,v,0) * cos(rv*pi)
    end do
    sher_0 = sher_0 - 1.D0
    sher_1 = sher_1 - 1.D0  
    write (30,*) 'sher_0 = ',sher_0
    write (30,*) 'sher_1 = ',sher_1 

    !Average Nusselt
    Nusselt_0=0.D0
    Nusselt_1=0.D0
    do s=1,ns
       rs = dble(s)
       Nusselt_0 = Nusselt_0 + pi * rs * GG(s,0,0)
       Nusselt_1 = Nusselt_1 + pi * rs * GG(s,0,0) * cos(rs*pi)
    end do
    Nusselt_0 = Nusselt_0 - 1.D0
    Nusselt_1 = Nusselt_1 - 1.D0  
    write (30,*) 'Nusselt_0 = ',Nusselt_0
    write (30,*) 'Nusselt_1 = ',Nusselt_1 
  end subroutine cal_Phi_C_vx_vy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end MODULE plot_module
