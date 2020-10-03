MODULE jacobian_module
  use com_module
  use residual_module
contains
  Subroutine Cal_derv_res(X,Jac,cjdaspk)
    use in_data
    implicit none 
    real*8 , dimension(:)    ,intent(in)  :: X(nve+nvg)
    real*8 , dimension(:,:)  ,intent(out) :: Jac(nve+nvg,nve+nvg)
    real*8                   , intent (in) :: cjdaspk
    integer                               :: cI,cJ,cK,cL,cM,cN,cU,cV,cW
    integer                               :: cS,cP,cT,i,j,k,l,m,n,u,v,w,s,&
         &p,t,sifx,sjfx,sify,sjfy,sitr,sjtr,siEn,sjEn
    real*8                                :: TermFy1,termFy2,termFy3,TermFx1,termFx2,termFx3   
    real*8 , dimension(:,:)               :: dAA_dEE(nva,nve),dBB_dEE(nvb,nve),dAA_dGG(nva,nvg),dBB_dGG(nvb,nvg)
    real*8 , dimension(:,:)               :: dFTR_dEE(nve,nve),dFTR_dAA(nve,nva),dFTR_dBB(nve,nvb),dFTR_dGG(nve,nvg)
    real*8 , dimension(:,:)               :: dFEn_dGG(nvg,nvg),dFEn_dBB(nvg,nvb),dFEn_dAA(nvg,nva)
    real*8 , dimension(:,:)               :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn),EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt)
    real*8                                :: TermFX4_Euvw,TermFX6_Gspt,TermFY4_Euvw,TermFy5_Gspt
    real*8                                :: TermTR1234_Aijk,TermTR9_Aijk,TermTR1234_Blmn,TermTR1234_Euvw,TermTR567_Euvw
    real*8                                :: TermEn1234_Aijk,TermEn9_Blmn,TermEn1234_Blmn,TermEn1234_Gspt,TermEn567_Gspt
    real*8 , dimension(:)                 :: der_yprim(nve+nvg)



    call vect_to_Matr (X,EE,GG)
    call Calc_AA(EE,GG,AA)
    call Calc_BB(EE,GG,BB)


    do cI=0,ni
       do cJ=1,nj
          do cK=1,nk
             sifx = CI * (nj*nk) + (Cj -1)* nk + CK  
             call cal_TermFX123 (CI,CJ,Ck,TermFx1,termFx2,termFx3)    

             do u=0,nu
                do v=1,nv
                   do w=0,nw
                      call der_TermFX4_Euvw (TermFX4_Euvw,cI,cJ,cK,u,v,w)
                      sjfx =  U * nv * (nw+1) + (V-1) * (nw+1) + W + 1
                      dAA_dEE (sifx,sjfx)= TermFX4_Euvw * Ngrav/ (TermFx1+TermFx2+TermFx3)
                   end do
                end do
             end do


             do s=1,ns
                do p=0,np
                   do t=0,nt
                      call der_TermFX6_Gspt (TermFX6_Gspt,cI,cJ,cK,s,p,t)
                      sjfx =  (S-1) * (np+1) * (nt+1) + P * (nt+1) + T + 1
                      dAA_dGG (sifx,sjfx)= (TermFX6_Gspt)/ (TermFx1+TermFx2+TermFx3)

                   end do
                end do
             end do

          end do
       end do
    end do

    ! dBB_dEE and dBB_dGG

    do cL=1,nl
       do cM=0,nm
          do cN=1,nn

             sify = (cL-1) * (nm+1) * nn + cM * nn + cN
             call cal_TermFY123 (cL,cM,CN,TermFy1,termFy2,termFy3)

             do u=0,nu
                do v=1,nv
                   do w=0,nw
                      call der_TermFy4_Euvw (TermFy4_Euvw,cL,cM,cN,u,v,w)
                      sjfy = U * nv * (nw+1) + (V-1) * (nw+1) + W + 1
                      dBB_dEE(sify,sjfy)= -TermFY4_Euvw*Ngrav/(TermFy1 + TermFy2 +TermFy3 )

                   end do
                end do
             end do


             do s=1,ns
                do p=0,np
                   do t=0,nt
                      call der_TermFy5_Gspt (TermFy5_Gspt,cL,cM,cN,s,p,t)
                      sjfy = (S-1) * (np+1) * (nt+1) + P * (nt+1) + T + 1
                      dBB_dGG (sify,sjfy)= -TermFY5_Gspt/(TermFy1 + TermFy2 +TermFy3 )

                   end do
                end do
             end do


          end do
       end do
    end do

    ! dFTR_dAA, dFTR_dBB, dFTR_dEE

    do cU=0,nu
       do cV=1,nv
          do cW=0,nw

             sitr = cU * nv * (nw+1) + (cV-1) * (nw+1) + cW + 1

             do i=0,ni 
                do j=1,nj 
                   do k=1,nk 
                      call der_TermTR1234_Aijk (TermTR1234_Aijk,EE,cU,cV,cW,i,j,k)

                      call der_TermTR9_Aijk (TermTR9_Aijk,cU,cV,cW,i,j,k)

                      sjtr = I * (nj*nk) + (J-1) * nk + K 
                      dFTR_dAA (sitr,sjtr) = TermTR1234_Aijk - TermTR9_Aijk  

                   end do
                end do
             end do

             do l=1,nl
                do m=0,nm
                   do n=1,nn
                      call der_TermTR1234_Blmn (TermTR1234_Blmn,EE,cU,cV,cW,l,m,n)

                      sjtr =  (L-1) * (nm+1) * nn + M * nn + N
                      dFTR_dBB (sitr,sjtr) = TermTR1234_Blmn 
                   end do
                end do
             end do

             do u=0,nu
                do v=1,nv
                   do w=0,nw
                      call der_TermTR1234_Euvw (TermTR1234_Euvw,AA,BB,cU,cV,cW,u,v,w)
                      call der_TermTR567_Euvw (TermTR567_Euvw,cU,cV,cW,u,v,w) 
                      sjtr =  U * nv * (nw+1) + (V-1) * (nw+1) + W + 1 
                      dFTR_dEE (sitr,sjtr) = TermTR1234_Euvw - TermTR567_Euvw/Lewis

                   end do
                end do
             end do

             call der_TermTR0_t_dEEdtuvw (cU,cV,cW,der_yprim (sitr))


          end do
       end do
    end do

    jac (1:nve,1:nve) =  dFTR_dEE + matmul(dFTR_dAA,dAA_dEE) + matmul(dFTR_dBB,dBB_dEE)
    jac (1:nve,nve+1:nve+nvg) = matmul(dFTR_dAA,dAA_dGG) + matmul(dFTR_dBB,dBB_dGG)


    do cS=1,ns
       do cp=0,np
          do cT=0,nt

             siEn = (cS-1) * (np+1) * (nt+1) + cP * (nt+1) + cT + 1   

             do i=0,ni 
                do j=1,nj 
                   do k=1,nk 
                      call der_TermEn1234_Aijk (TermEn1234_Aijk,GG,cS,cP,cT,i,j,k)

                      sjEn = I * (nj*nk) + (J-1) * nk + K 
                      dFEn_dAA (siEn,sjEn) = TermEn1234_Aijk   

                   end do
                end do
             end do


             do l=1,nl
                do m=0,nm
                   do n=1,nn
                      call der_TermEn1234_Blmn (TermEn1234_Blmn,GG,cS,cP,cT,l,m,n)
                      call der_TermEn9_Blmn (TermEn9_Blmn,cS,cP,cT,l,m,n)
                      sjEn = (L-1) * (nm+1) * nn + M * nn + N
                      dFEn_dBB (siEn,sjEn) =  TermEn1234_Blmn + TermEn9_Blmn 
                   end do
                end do
             end do


             do s=1,ns
                do p=0,np
                   do t=0,nt
                      call der_TermEn1234_Gspt (TermEn1234_Gspt,AA,BB,cS,cP,cT,s,p,t)
                      call der_TermEn567_Gspt (TermEn567_Gspt,cS,cP,cT,s,p,t) 
                      sjEn = (S-1) * (np+1) * (nt+1) + P * (nt+1) + T + 1 
                      dFEn_dGG (siEn,sjEn)= TermEn1234_Gspt - TermEn567_Gspt 

                   end do
                end do
             end do

             call der_TermEn0_t_dGGdtspt (cS,cP,cT,der_yprim (siEn+nve))

          end do
       end do
    end do


    jac (nve+1:nve+nvg,1:nve)         = matmul(dFEn_dAA,dAA_dEE) + matmul(dFEn_dBB,dBB_dEE)
    jac (nve+1:nve+nvg,nve+1:nve+nvg) = dFEn_dGG  + matmul(dFEn_dAA,dAA_dGG) + matmul(dFEn_dBB,dBB_dGG)

    do i = 1, nve + nvg
       jac(i,i) = jac(i,i)+ der_yprim (i) * cjdaspk
    end do


    return 
  end Subroutine Cal_derv_res
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermFX123_Aijk (TermFX123_Aijk,cI,cJ,cK,i,j,k)
    use in_data  
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK,i,j,k
    real*8                       , intent (out) :: termFx123_Aijk
    real*8                                      :: cIr,cJr,cKr,TermFX1_Aijk,TermFX2_Aijk,TermFX3_Aijk

    cIr = dble(cI)
    cJr = dble(cJ)
    cKr = dble(cK)

    TermFX123_Aijk=0.d0
    TermFX1_Aijk=0.d0
    TermFX2_Aijk=0.d0
    TermFX3_Aijk=0.d0

    if ((cI==i).and.(cJ==j).and.(cK==k)) then
       termFX1_Aijk = -(pi**2.d0) *  cIr**(2.d0)
       if (cI.eq.0) then  
          termFX2_Aijk = -2.d0* (pi**2.d0) *cJr**(2.d0)
          termFX3_Aijk = -2.d0* (pi**2.d0) *cKr**(2.d0)
       else 
          termFX2_Aijk = -(pi**2.d0) *cJr**(2.d0)
          termFX3_Aijk = -(pi**2.d0) *cKr**(2.d0)
       end if
       termFX123_Aijk = termFX1_Aijk + termFX2_Aijk + termFX3_Aijk   
    else
       termFX123_Aijk=0.d0      
    end if


  end subroutine der_TermFX123_Aijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermFX4_Euvw (TermFX4_Euvw,cI,cJ,cK,u,v,w)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK,u,v,w
    real*8                       , intent (out) :: termFx4_Euvw
    real*8                                      :: vr

    termFx4_Euvw = 0.d0 

    if (cI==u)then
       vr = dble(v)
       termFx4_Euvw = vr * matphi(cJ,v) * matphi(cK,w)  

       if (cI.eq.0) then 
          termFx4_Euvw = 2.d0*termFx4_Euvw * Ra / pi
       else 
          termFx4_Euvw = termFx4_Euvw * Ra /  pi
       end if
    else
       termFx4_Euvw=0.d0   
    end if

  end subroutine der_TermFX4_Euvw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermFX6_Gspt (TermFX6_Gspt,cI,cJ,cK,s,p,t)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK,s,p,t
    real*8                       , intent (out) :: termFx6_Gspt
    real*8                                      :: pr


    termFx6_Gspt = 0.d0 
    if (cJ .eq. p) then
       termFx6_Gspt=termFx6_Gspt-dble(p)*matphi(s,cI)*matphi(cK,t)
    end if




    TermFX6_Gspt = TermFX6_Gspt * Ra /  pi
  end subroutine der_TermFX6_Gspt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermFY123_Blmn (TermFY123_Blmn,cL,cM,cN,l,m,n)
    use in_data  
    implicit none 
    integer                      , intent (in)  :: cL,cM,cN,l,m,n
    real*8                       , intent (out) :: TermFY123_Blmn
    real*8                                      :: cLr,cMr,cNr,termFy1_Blmn,termFy2_Blmn,termFy3_Blmn


    cLr = dble(cL)
    cMr = dble(cM)
    cNr = dble(cN)

    TermFY123_Blmn=0.d0
    TermFY1_Blmn=0.d0
    TermFY2_Blmn=0.d0
    TermFY3_Blmn=0.d0

    if ((cL==l).and.(cM==m).and.(cN==n)) then
       termFy2_Blmn = - (pi*cMr)**2.d0

       if (cM.eq.0) then  
          termFy1_Blmn = -2.d0* (pi**2.d0) *  cLr**(2.d0)
          termFy3_Blmn = -2.d0* (pi**2.d0) *  cNr**(2.d0)
       else 
          termFy1_Blmn = - (pi**2.d0) *  cLr**(2.d0)
          termFy3_Blmn = - (pi**2.d0) *  cNr**(2.d0)
       end if
       termFY123_Blmn = termFY1_Blmn + termFY2_Blmn + termFY3_Blmn
    else
       TermFY123_Blmn = 0.D0
    end if


  end subroutine der_TermFY123_Blmn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermFy4_Euvw (TermFy4_Euvw,cL,cM,cN,u,v,w)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cL,cM,cN,u,v,w
    real*8                       , intent (out) :: TermFy4_Euvw
    real*8                                      :: cLr

    termFy4_Euvw = 0.d0 
    cLr      = dble (cL)

    if (cL==u)then
       termFy4_Euvw = matphi(v,cM) * matphi(cN,w)
    else 
       termFy4_Euvw = 0.d0 
    end if
    termFy4_Euvw = -cLr * Ra * termFy4_Euvw  / pi

  end subroutine der_TermFy4_Euvw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermFy5_Gspt(TermFy5_Gspt,cL,cM,cN,s,p,t)
    use in_data
    use matrices
    implicit none
    integer                      , intent (in)  :: cL,cM,cN,s,p,t
    real*8                       , intent (out) :: termFy5_Gspt


    termFy5_Gspt = 0.d0 

    if (cM .eq. p) then 

       termFy5_Gspt = termFy5_Gspt + dble(s)* matphi(cL,s) * matphi(cN,t)
       if (cM.eq.0) then 
          termFy5_Gspt = 2.d0*termFy5_Gspt * Ra / pi
       else 
          termFy5_Gspt = termFy5_Gspt * Ra /  pi
       end if
    else 
       termFy5_Gspt = 0.d0
    end if



  end subroutine der_TermFy5_Gspt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermTR1234_Aijk (TermTREn1234_Aijk,EE,cU,cV,cW,i,j,k)
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW,i,j,k
    real*8    , dimension (:,:)  , intent (in)  :: EE(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: TermTREn1234_Aijk
    real*8                                      :: TermTREn2_Aijk,TermTREn4_Aijk
    INTEGER :: u,v,w

    TermTREn1234_Aijk=0.D0
    termTREn2_Aijk=0.d0
    termTREn4_Aijk=0.d0

    if (cU .eq. 0)then
       u=cU-i
       if ((u .le. nu) .and. (u .ge. 0))then

          v=cV-j
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif

          v=j-cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
          v=j+cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
       end if
       u=cU+i
       if ((u .le. nu .and. u .ge. 0))then
          v=cV-j
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif

          v=j-cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
          v=j+cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/4.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/4.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
       end if
    else
       u=cU-i
       if ((u .le. nu) .and. (u .ge. 0))then

          v=cV-j
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif

          v=j-cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
          v=j+cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
       end if
       u=i-cU
       if (( u .ge. 0) .and. ( u.le. nu)) then
          v=cV-j
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif

          v=j-cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
          v=j+cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
       end if
       u=cU+i
       if ((u .ge. 0 ) .and. (u.le. nu)) then
          v=cV-j
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif

          v=j-cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk+(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
          v=j+cV
          if((v.le.nv).and.(v.ge.1))then
             w=cW-k
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k+cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
             w=k-cW
             if((w.le.nw).and.(w.ge.0))then
                termTREn2_Aijk=termTREn2_Aijk-(pi**2.d0)/8.0d0*dble(v)*dble(k)*EE(u,v,w)
                termTREn4_Aijk=termTREn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(w)*EE(u,v,w)
             endif
          endif
       end if
    end if

    TermTREn1234_Aijk= - termTREn4_Aijk + termTREn2_Aijk

  end subroutine der_TermTR1234_Aijk



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermTR1234_Blmn (TermTREn1234_Blmn,EE,cU,cV,cW,l,m,n)
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW,l,m,n
    real*8    , dimension (:,:)  , intent (in)  :: EE(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: TermTREn1234_Blmn
    real*8                                      :: termTREn1_Blmn,termTREn3_Blmn
    INTEGER :: u,v,w

    TermTREn1234_Blmn=0.D0
    termTREn1_Blmn=0.d0
    termTREn3_Blmn=0.d0


    if ( cU .eq. 0) then
       u=  cU - l
       if ( (u .ge. 0) .and. (u .le. nu)) then
          v= cV-m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=cV+m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=m-cV
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
       end if
       u= cU + l
       if (( u .ge. 0) .and. ( u.le. nu)) then
          v= cV-m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=cV+m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=m-cV
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
       end if
    else 
       u= cU-l
       if ( (u .ge. 0) .and. (u .le. nu)) then
          v= cV-m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=cV+m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=m-cV
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
       end if
       u= cU + l
       if (( u .ge. 0) .and. ( u.le. nu)) then
          v= cV-m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=cV+m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=m-cV
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
       end if
       u=l-cU
       if (( u .ge. 0) .and. ( u.le. nu)) then
          v= cV-m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=cV+m
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
          v=m-cV
          if ( (v .ge. 1) .and. (v .le. nv)) then
             w= cW- n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=cW+n
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
             w=n-cW  
             if ( (w .ge. 0) .and. (w .le. nw)) then
                termTREn1_Blmn=termTREn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(u)*EE(u,v,w)
                termTREn3_Blmn= termTREn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(w)*EE(u,v,w)
             end if
          end if
       end if
    end if


    TermTREn1234_Blmn = TermTREn3_Blmn -TermTREn1_Blmn
  end subroutine der_TermTR1234_Blmn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermTR1234_Euvw (TermTREn1234_EGuvw,AA,BB,cU,cV,cW,u,v,w)
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW,u,v,w
    real*8    , dimension (:,:)  , intent (in)  :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn)
    real*8                       , intent (out) :: TermTREn1234_EGuvw
    real*8                                      :: termTREn1_EGuvw,termTREn2_EGuvw,termTREn3_EGuvw,termTREn4_EGuvw
    INTEGER :: i, j, k, l, m, n

    TermTREn1234_EGuvw=0.D0
    termTREn1_EGuvw=0.d0
    termTREn2_EGuvw=0.d0
    termTREn3_EGuvw=0.d0
    termTREn4_EGuvw=0.d0
    if (cU .eq. 0) then
       i=cU-u
       if ((i.le. ni) .and. (i.ge.0)) then
          j=cV-v  !!!d5
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= cV+v  !!!d6
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= v-cV  !!!-d7
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
       end if
       i=u-cU
       if ((i.le. ni) .and. (i.ge.0)) then
          j=cV-v  !!!d5
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= cV+v  !!!d6
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= v-cV  !!!-d7
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/4.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/4.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
       end if
    else
       i=cU-u
       if ((i.le. ni) .and. (i.ge.0)) then
          j=cV-v  !!!d5
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= cV+v  !!!d6
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= v-cV  !!!-d7
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
       end if
       i=u-cU
       if ((i.le. ni) .and. (i.ge.0)) then
          j=cV-v  !!!d5
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= cV+v  !!!d6
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= v-cV  !!!-d7
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
       end if
       i=cU+u
       if ((i.le. ni) .and. (i.ge.0)) then
          j=cV-v  !!!d5
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= cV+v  !!!d6
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw+(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
          j= v-cV  !!!-d7
          if ((j.le. nj) .and. (j .ge.1)) then
             k=cW-w
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw+(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w-cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
             k=w+cW
             if ((k .le. nk) .and. (k.ge.1)) then
                termTREn2_EGuvw=termTREn2_EGuvw-(pi**2.d0)/8.0d0*dble(k)*dble(v)*AA(i,j,k)
                termTREn4_EGuvw=termTREn4_EGuvw-(pi**2.d0)/8.0d0*dble(j)*dble(w)*AA(i,j,k)
             end if
          end if
       end if
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!end of the calculation of the derivative term2,4 with respect to E
    if (cU .eq. 0) then
       l=cU-u
       if((l .le.nl) .and. (l .ge.1)) then
          m=cV-v
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v-cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v+cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
       end if
       l=u-cU
       if((l .le.nl) .and. (l .ge.1)) then
          m=cV-v
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v-cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v+cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/4.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/4.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
       end if
    else
       l=cU-u
       if((l .le.nl) .and. (l .ge.1)) then
          m=cV-v
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v-cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v+cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
       end if
       l=u-cU
       if((l .le.nl) .and. (l .ge.1)) then
          m=cV-v
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v-cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v+cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
       end if
       l=u+cU
       if((l .le.nl) .and. (l .ge.1)) then
          m=cV-v
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v-cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw-(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
          m=v+cV
          if ((m .le. nm) .and. (m .ge.0)) then
             n=cW-w
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw-(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w-cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
             n=w+cW
             if ((n .le.nn) .and. (n .ge. 1)) then
                termTREn1_EGuvw=termTREn1_EGuvw+(pi**2.d0)/8.0d0*dble(n)*dble(u)*BB(l,m,n)
                termTREn3_EGuvw=termTREn3_EGuvw+(pi**2.d0)/8.0d0*dble(w)*dble(l)*BB(l,m,n)
             end if
          end if
       end if
    end if


    TermTREn1234_EGuvw= - TermTREn1_EGuvw + TermTREn2_EGuvw + TermTREn3_EGuvw - TermTREn4_EGuvw   

  end subroutine der_TermTR1234_Euvw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermTR567_Euvw (TermTR567_Euvw,cU,cV,cW,u,v,w)       
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW,u,v,w
    real*8                       , intent (out) :: TermTR567_Euvw
    real*8                                      :: cUr,cVr,cWr,termTR5_Euvw,termTR6_Euvw,termTR7_Euvw

    TermTR567_Euvw=0.D0
    TermTR5_Euvw=0.D0
    TermTR6_Euvw=0.D0
    TermTR7_Euvw=0.D0

    cUr=dble(cU)
    cVr=dble(cV)
    cWr=dble(cW)

    if ((u==cU).and.(v==cV).and.(w==cW)) then
       if (cW.eq.0) then
          termTR5_Euvw=-2.d0*(pi**2.d0)*(cUr**2.d0)
       else
          termTR5_Euvw=-(pi**2.d0)*(cUr**2.d0) 
       endif

       if ((cW.eq.0).and.(cU.eq.0)) then
          termTR6_Euvw=-4.d0*(pi**2.d0)*(cVr**2.d0) 
       elseif ((cW.ne.0).and.(cU.ne.0)) then
          termTR6_Euvw=-(pi**2.d0)*(cVr**2.d0)  
       else
          termTR6_Euvw=-2.d0*(pi**2.d0)*(cVr**2.d0)
       endif

       if (cU.eq.0) then
          termTR7_Euvw=-2.d0*(pi**2.d0)*(cWr**2.d0)
       else
          termTR7_Euvw=-(pi**2.d0)*(cWr**2.d0)
       endif
    endif

    TermTR567_Euvw = TermTR5_Euvw + TermTR6_Euvw + TermTR7_Euvw 

  end subroutine der_TermTR567_Euvw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermTR9_Aijk (TermTR9_Aijk,cU,cV,cW,i,j,k)       
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW,i,j,k
    real*8                       , intent (out) :: termTR9_Aijk

    real*8 :: d1, d2, d3
    termTR9_Aijk=0.D0
    call delta_funct(i,cU,d1)
    call delta_funct(j,cV,d2)
    call delta_funct(k,cW,d3)
    if (cU .eq. 0) then
       termTR9_Aijk=termTR9_Aijk+2*pi*dble(k)* d1*d2*d3
    else 
       termTR9_Aijk=termTR9_Aijk+pi*dble(k)* d1*d2*d3
    end if


  end subroutine der_TermTR9_Aijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermEn1234_Aijk (TermEn1234_Aijk,GG,cS,cP,cT,i,j,k)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT,i,j,k
    real*8    , dimension (:,:)  , intent (in)  :: GG(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: TermEn1234_Aijk
    real*8                                      :: TermEn2_Aijk,TermEn4_Aijk
    INTEGER :: s,p,t

    TermEn1234_Aijk=0.D0
    termEn2_Aijk=0.d0
    termEn4_Aijk=0.d0
    s=cS-i
    if ((s .le.ns) .and. (s .ge. 1)) then
       p= j-cP
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
       p=cP+j
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
       p=cP-j  
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if

    end if
    s=cS+i
    if ((s .le.ns) .and. (s .ge. 1)) then
       p= j-cP
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
       p=cP+j
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
       p=cP-j  
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
    end if
    s=i-cS
    if ((s .le.ns) .and. (s .ge. 1)) then
       p= j-cP
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
       p=cP+j
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk+(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if
       p=cP-j  
       if ((p .le. np) .and. (p .ge. 0)) then
          t= cT-k
          if ((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk-(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t= k-cT
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
          t=cT+k 
          if((t .le. nt) .and. (t .ge. 0)) then
             termEn2_Aijk=termEn2_Aijk-(pi**2.d0)/8.0d0*dble(k)*dble(p)*GG(s,p,t)
             termEn4_Aijk=termEn4_Aijk+(pi**2.d0)/8.0d0*dble(j)*dble(t)*GG(s,p,t)
          end if
       end if

    end if
    TermEn1234_Aijk=termEn2_Aijk-termEn4_Aijk
  end subroutine der_TermEn1234_Aijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermEn1234_Blmn (TermEn1234_Blmn,GG,cS,cP,cT,l,m,n)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT,l,m,n
    real*8    , dimension (:,:)  , intent (in)  :: GG(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: TermEn1234_Blmn
    real*8                                      :: TermEn1_Blmn,TermEn3_Blmn
    INTEGER :: s,p,t

    TermEn1234_Blmn=0.D0
    TermEn1_Blmn=0.d0
    TermEn3_Blmn=0.d0
    if (cP .eq. 0) then
       s=cS-l
       if ((s .le.ns) .and. (s.ge.1)) then
          p=cP-m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=cP+m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
       end if
       s=l-cS  
       if ((s .le.ns) .and. (s.ge.1)) then
          p=cP-m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=cP+m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
       end if
       s=l+cS
       if ((s .le.ns) .and. (s.ge.1)) then
          p=cP-m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=cP+m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/4.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/4.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
       end if
    else
       s=cS-l
       if ((s .le.ns) .and. (s.ge.1)) then
          p=cP-m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=cP+m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=m-cP
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
       end if
       s=l-cS
       if ((s .le.ns) .and. (s.ge.1)) then
          p=cP-m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=cP+m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=m-cP
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn+(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
       end if
       s=l+cS
       if ((s .le.ns) .and. (s.ge.1)) then
          p=cP-m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=cP+m
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
          p=m-cP
          if ((p .le. np) .and. (p .ge.0)) then
             t=cT-n
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn+(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n-cT
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
             t=n+cT  
             if ((t .le. nt) .and. (t.ge.0)) then
                TermEn1_Blmn=TermEn1_Blmn-(pi**2.d0)/8.0d0*dble(n)*dble(s)*GG(s,p,t)
                TermEn3_Blmn=TermEn3_Blmn-(pi**2.d0)/8.0d0*dble(l)*dble(t)*GG(s,p,t)
             end if
          end if
       end if
    end if
    TermEn1234_Blmn=TermEn3_Blmn-TermEn1_Blmn 
  end subroutine der_TermEn1234_Blmn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermEn9_Blmn (TermEn9_Blmn,cS,cP,cT,l,m,n)   
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT,l,m,n
    real*8                       , intent (out) :: termEn9_Blmn
    real*8 :: d1, d2, d3
    termEn9_Blmn=0.D0
    call delta_funct(l,cS,d1)
    call delta_funct(m,cP,d2)
    call delta_funct(n,cT,d3)
    if (cP .eq. 0) then
       termEn9_Blmn=termEn9_Blmn+2* dble(n) * d1*d2*d3
    else 
       termEn9_Blmn=termEn9_Blmn+  dble(n) * d1*d2*d3
    end if
    termEn9_Blmn=termEn9_Blmn*pi

  end subroutine der_TermEn9_Blmn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermEn1234_Gspt (TermEn1234_Gspt,AA,BB,cS,cP,cT,s,p,t)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT,s,p,t
    real*8    , dimension (:,:)  , intent (in)  :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn)
    real*8                       , intent (out) :: TermEn1234_Gspt
    real*8                                      :: termEn1_Gspt,termEn2_Gspt,termEn3_Gspt,termEn4_Gspt
    INTEGER :: i, j, k, l, m, n

    TermEn1234_Gspt=0.D0
    TermEn1_Gspt=0.d0
    TermEn2_Gspt=0.d0
    TermEn3_Gspt=0.d0
    TermEn4_Gspt=0.d0
    i=cS-s
    if((i .le. ni) .and. (i .ge. 0)) then
       j=cP+p
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
       j=p-cP
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
       j=cP-p
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
    end if
    i=s-cS
    if((i .le. ni) .and. (i .ge. 0)) then
       j=cP+p
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
       j=p-cP
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
       j=cP-p
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
    end if
    i=s+cS  
    if((i .le. ni) .and. (i .ge. 0)) then
       j=cP+p
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
       j=p-cP
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt+(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
       j=cP-p
       if ((j .le. nj) .and. (j .ge. 1)) then
          k=cT-t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt-(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=cT+t
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
          k=t-cT
          if((k .le. nk) .and. (k .ge. 1)) then
             TermEn2_Gspt=TermEn2_Gspt-(pi**2.d0)/8.0d0*dble(k)*dble(p)*AA(i,j,k)
             TermEn4_Gspt=TermEn4_Gspt+(pi**2.d0)/8.0d0*dble(j)*dble(t)*AA(i,j,k)
          end if
       end if
    end if
    if( cP .eq. 0) then
       l=cS-s
       if ((l .le. nl) .and. (l.ge. 1)) then
          m=cP-p
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p-cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
       end if
       l=cS+s
       if ((l .le. nl) .and. (l.ge. 1)) then
          m=cP-p
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p-cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
       end if
       l=s-cS
       if ((l .le. nl) .and. (l.ge. 1)) then
          m=cP-p
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p-cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/4.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/4.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
       end if
    else
       l=cS-s
       if ((l .le. nl) .and. (l.ge. 1)) then
          m=cP-p
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p-cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p+cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
       end if
       l=cS+s
       if ((l .le. nl) .and. (l.ge. 1)) then
          m=cP-p
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p-cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p+cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt+(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
       end if
       l=s-cS 
       if ((l .le. nl) .and. (l.ge. 1)) then
          m=cP-p
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p-cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
          m=p+cP
          if ((m .le.nm) .and. (m .ge. 0)) then
             n=cT-t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt+(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=cT+t
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
             n=t-cT
             if ((n.le. nn) .and. (n .ge. 1)) then
                TermEn1_Gspt=TermEn1_Gspt-(pi**2.d0)/8.0d0*dble(n)*dble(s)*BB(l,m,n) 
                TermEn3_Gspt=TermEn3_Gspt-(pi**2.d0)/8.0d0*dble(l)*dble(t)*BB(l,m,n)  
             end if
          end if
       end if
    end if
    TermEn1234_Gspt= - TermEn1_Gspt + TermEn2_Gspt + TermEn3_Gspt - TermEn4_Gspt
  end subroutine der_TermEn1234_Gspt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine der_TermEn567_Gspt (TermEn567_Gspt,cS,cP,cT,s,p,t) 
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT,s,p,t
    real*8                       , intent (out) :: TermEn567_Gspt
    real*8                                      :: cPr,cSr,cTr,termEn5_Gspt,termEn6_Gspt,termEn7_Gspt

    TermEn567_Gspt=0.D0
    TermEn5_Gspt=0.D0
    TermEn6_Gspt=0.D0
    TermEn7_Gspt=0.D0

    cPr=dble(cP)
    cSr=dble(cS)
    cTr=dble(cT)

    if ((s==cS).and.(p==cP).and.(t==cT)) then
       if ((cP.eq.0).and.(cT.eq.0)) then
          termEn5_Gspt=-4.d0*(pi**2.d0)*(cSr**2.d0) 
       elseif ((cT.ne.0).and.(cP.ne.0)) then
          termEn5_Gspt=-(pi**2.d0)*(cSr**2.d0)  
       else
          termEn5_Gspt=-2.d0*(pi**2.d0)*(cSr**2.d0)
       endif

       if (cT.eq.0) then
          termEn6_Gspt=-2.d0*(pi**2.d0)*(cPr**2.d0)
       else
          termEn6_Gspt=-(pi**2.d0)*(cPr**2.d0) 
       endif

       if (cP.eq.0) then
          termEn7_Gspt=-2.d0*(pi**2.d0)*(cTr**2.d0)
       else
          termEn7_Gspt=-(pi**2.d0)*(cTr**2.d0)
       endif
    endif
    !
    TermEn567_Gspt = TermEn5_Gspt + TermEn6_Gspt + TermEn7_Gspt 

  end subroutine der_TermEn567_Gspt



  subroutine der_TermTR0_t_dEEdtuvw (cU,cV,cW,TermTR0_t_dEEdtuvw)
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8                       , intent (out) :: TermTR0_t_dEEdtuvw

    TermTR0_t_dEEdtuvw = poro 

    if (Cu.eq. 0 .and. CW .eq. 0) TermTR0_t_dEEdtuvw = 4.d0 * poro 
    if (Cu.eq. 0 .and. CW .ne. 0) TermTR0_t_dEEdtuvw = 2.d0 * poro 
    if (Cu.ne. 0 .and. CW .eq. 0) TermTR0_t_dEEdtuvw = 2.d0 * poro 

  end subroutine der_TermTR0_t_dEEdtuvw



!!!!!!!!!!!!!!!
  subroutine  der_TermEn0_t_dGGdtspt (cS,cP,cT,TermEn0_t_dGGdtspt)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT
    real*8                       , intent (out) :: TermEn0_t_dGGdtspt


    TermEn0_t_dGGdtspt = sigma 

    if (Cp.eq. 0 .and. CT .eq. 0) TermEn0_t_dGGdtspt = 4.d0 *  sigma 
    if (Cp.eq. 0 .and. CT .ne. 0) TermEn0_t_dGGdtspt = 2.d0 *  sigma 
    if (Cp.ne. 0 .and. CT .eq. 0) TermEn0_t_dGGdtspt = 2.d0 *  sigma 

  end subroutine der_TermEn0_t_dGGdtspt

end MODULE jacobian_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
