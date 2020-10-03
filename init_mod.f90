MODULE init_module
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_data_in_cal_coef()
    use in_data

    read(17,*) Ra
    read(17,*) lewis
    read(17,*) Ngrav
    read(17,*) poro , sigma
    read(17,*) ni
    read(17,*) nj
    read(17,*) nk
    read(17,*) nl
    read(17,*) nm
    read(17,*) nn
    read(17,*) nu
    read(17,*) nv
    read(17,*) nw
    read(17,*) ns
    read(17,*) np
    read(17,*) nt
    nva = (ni+1)*nj*nk
    nvb =  nl*(nm+1)*nn
    nve = (nu+1)*nv*(nw+1)
    nvg= ns*(np+1)*(nt+1)

    ntrucmax = ni
    if (nj.gt.ntrucmax) ntrucmax = nj 
    if (nk.gt.ntrucmax) ntrucmax = nk
    if (nl.gt.ntrucmax) ntrucmax = nl
    if (nm.gt.ntrucmax) ntrucmax = nm
    if (nn.gt.ntrucmax) ntrucmax = nn
    if (nu.gt.ntrucmax) ntrucmax = nu
    if (nv.gt.ntrucmax) ntrucmax = nv
    if (nw.gt.ntrucmax) ntrucmax = nw
    if (ns.gt.ntrucmax) ntrucmax = ns
    if (np.gt.ntrucmax) ntrucmax = np
    if (nt.gt.ntrucmax) ntrucmax = nt
    Read(49,*) ntime
    allocate  (Tout(ntime))
    Read(49,*) tout(1:ntime)


    RETURN
  end subroutine read_data_in_cal_coef
  !**************************************************************************************
  subroutine Init_Guess(EE,GG)
    use in_data  
    real*8 , dimension(:,:,:),intent(out)  :: EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt)
    integer                                :: co1,co2,co3,test

    print*,'Initial values from file (init_coef) ?    ***  1=yes / 0 = No'

    read(*,*) test

    if (test.eq.0) then 


       do co1=0,nu
          do co2=1,nv 
             do co3=0,nw
                EE(co1,co2,co3) =  0.2d0*co1-0.1d0*co2-0.2d0*co3
                !print*,co1,co2,co3,EE(co1,co2,co3)

             end do
          end do
       end do

       do co1=1,ns
          do co2=0,np 
             do co3=0,nt
                GG(co1,co2,co3) =  0.4d0*co1+0.1d0*co2+0.2d0*co3

             end do
          end do
       end do

    else 

       do co1=0,nu
          do co2=1,nv 
             do co3=0,nw
                read(40,*) EE(co1,co2,co3) 
             end do
          end do
       end do

       do co1=1,ns
          do co2=0,np 
             do co3=0,nt
                read(40,*) GG(co1,co2,co3) 
             end do
          end do
       end do


    end if

    RETURN

  end subroutine Init_Guess
  !**************************************************************************************
end MODULE init_module
