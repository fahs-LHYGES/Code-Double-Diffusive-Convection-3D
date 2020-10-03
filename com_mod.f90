MODULE com_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  subroutine matr_to_vect(EE,GG,vectABE)
    use in_data  
    implicit none 
    real*8 , dimension(:,:),intent(in)  :: EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt)
    real*8 , dimension(:)  ,intent(out) :: vectABE(nve+nvg)
    integer                             :: co1,co2,co3,s 

    s=1

    do co1 = 0,nu
       do co2= 1,nv
          do co3= 0,nw
             vectABE (s) = EE(co1,co2,co3)
             s = s + 1 
          end do
       end do
    end do

    do co1=1,ns
       do co2=0,np 
          do co3=0,nt
             vectABE (s) = GG(co1,co2,co3)
             s = s + 1 
          end do
       end do
    end do


    RETURN
  end subroutine matr_to_vect
  !**************************************************************************************
  subroutine vect_to_Matr(vectABE,EE,GG)
    use in_data
    implicit none 
    real*8 , dimension(:,:),intent(out)  :: EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt)
    real*8 , dimension(:)  ,intent(in)   :: vectABE(nve+nvg)

    integer :: co1,co2,co3,s 

    s=1

    do co1 = 0,nu
       do co2= 1,nv
          do co3= 0,nw
             EE(co1,co2,co3) = vectABE (s)
             s = s + 1
          end do
       end do
    end do

    do co1=1,ns
       do co2=0,np 
          do co3=0,nt
             GG(co1,co2,co3) = vectABE (s)
             s = s + 1
          end do
       end do
    end do


    RETURN
  end subroutine vect_to_Matr
  !**************************************************************************************
  subroutine cal_matphi ()
    use matrices
    use in_data 
    integer :: co1,co2
    real*8  :: rco1,rco2,t1,t2,t3,t4,t5,t6,t7,t8

    do co1=0,ntrucmax 
       do co2 = 0,ntrucmax 
          if (co1 .eq. co2) then 
             matphi(co1,co2) = 0.d0
          else   
             rco1 = dble (co1) 
             rco2 = dble (co2)
             t1 = rco1 + rco2
             t2 = (-1.d0)**(t1)
             t3 = 1.d0 - t2      
             t4 = t3/t1

             t5 = rco1 - rco2
             t6 = (-1.d0)**t5
             t7 = 1.d0 - t6
             t8 = t7/t5

             matphi(co1,co2) = t4 + t8
          end if
       end do
    end do

    return 
  end subroutine cal_matphi
  !**************************************************************************************
  subroutine delta_funct(i,j,dij)
    implicit none 
    integer, intent(in) :: i,j 
    real*8, intent(out) :: dij


    dij = 0.d0 
    if (i.eq.j) then 
       dij= 1.d0
    end if

  end subroutine delta_funct
end MODULE com_module
