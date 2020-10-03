MODULE residual_module
  use com_module
contains
  !**************************************************************************************
  Subroutine Cal_res(dxdt,X,F)
    use in_data  
    use OMP_LIB
    implicit none 
    real*8 , dimension(:)  ,intent(in)  :: X(nve+nvg),dXdt(nve+nvg)
    real*8 , dimension(:)  ,intent(out) :: F(nve+nvg)
    integer                             :: sfx,sfy,sEn,str,cI,cJ,ck,cL,cM,cN,cU,cV,cW,cS,cP,cT,i
    real*8                              :: norm,TermFx1,TermFx2,TermFx3,TermFx4,TermFx5,TermFx6
    real*8                              :: TermFy1,TermFy2,TermFy3,TermFy4,TermFy5,TermFy6
    real*8                              :: TermTR1,TermTR2,TermTR3,TermTR4,TermTR5,TermTR6,TermTR7,TermTR8,TermTR9,TermTR10
    real*8                              :: TermE1,TermE2,TermE3
    real*8                              :: TermE4,TermE5,TermE6
    real*8                              :: TermE7,TermE9,termTR0_t,TermEn0_t
    real*8 , dimension(:,:)             :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn)
    real*8 , dimension(:,:)             :: EE(0:nu,1:nv,0:nw),GG(1:ns,0:np,0:nt)
    real*8 , dimension(:,:)             :: dEEdt(0:nu,1:nv,0:nw),dGGdt(1:ns,0:np,0:nt)
    integer                             :: i_min,i_max,rang,nb_taches,l_min,l_max,u_min,u_max,s_min,s_max

    norm  = 0.d0
    F(:) = 0.d0 
    compt = compt +  1 

    call vect_to_Matr (X,EE,GG)
    call vect_to_Matr (dXdt,dEEdt,dGGdt)

    !$OMP PARALLEL PRIVATE(TermFx1,termFx2,termFx3,termFx4,termFx5,sfx,i_min,i_max,rang,TermFy1,termFy2,termFy3,termFy4,termFy5,sfy)
    !$OMP PARALLEL PRIVATE(l_min,l_max,termTR1,termTR2,termTR3,termTR4,termTR5,termTR6,termTR7,termTR8,termTR9,termTR10,str,u_min)
    !$OMP PARALLEL PRIVATE(u_max)
    rang= OMP_GET_THREAD_NUM () ; nb_taches= OMP_GET_NUM_THREADS () ; i_min=ni ; i_max=0;
    !$OMP DO SCHEDULE(STATIC,(ni+1)/nb_taches)

    do cI=0,ni
       i_min=min(i_min,CI) ; i_max=max(i_max,CI)
       do cJ=1,nj
          do cK=1,nk
             call cal_TermFX123 (cI,cJ,Ck,TermFx1,termFx2,termFx3) 
             call cal_TermFX4 (EE,cI,cJ,cK,termFx4)
             call cal_TermFX5 (cI,cJ,cK,termFx5)
             call cal_TermFX6 (GG,cI,cJ,cK,termFx6)
             AA (CI,Cj,Ck ) = (Ngrav*TermFx4 -Ngrav*TermFx5 +TermFx6 )/ (TermFx1+TermFx2+TermFx3)  
             ! print*,'AA[',ci,',',cj,',',ck,']:=',AA (CI,Cj,Ck ),';'
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    l_min=nl ; l_max=1;
    !$OMP DO SCHEDULE(STATIC,nl/nb_taches)
    do cL=1,nl
       l_min=min(l_min,CL) ; l_max=max(l_max,CL)
       do cM=0,nm
          do cN=1,nn
             call cal_TermFY123 (cL,cM,CN,TermFy1,termFy2,termFy3)
             call cal_TermFy4 (EE,cL,cM,cN,termFy4)
             call cal_TermFy5 (GG,cL,cM,cN,termFy5)
             call cal_TermFy6 (cL,cM,cN,termFy6)
             sfy = nva + (cL-1) * (nm+1) * nn + cM * nn + cN
             BB(Cl,CM,CN) = (-Ngrav * TermFy4 - TermFy5 +TermFy6)/(TermFy1 + TermFy2 +TermFy3 )
             ! print*,'BB[',cl,',',cm,',',cn,']:=',BB (Cl,Cm,Cn ),';'
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    u_min=nu ; u_max=0;
    !$OMP DO SCHEDULE(STATIC,(nu+1)/nb_taches)
    do cU=0,nu
       u_min=min(u_min,CU) ; u_max=max(u_max,CU)
       do cV=1,nv
          do cW=0,nw
             call cal_TermTR0_t (dEEdt,cU,cV,cW,termTR0_t)
             call cal_TermTR1234 (AA,BB,EE,cU,cV,cW,termTR1,termTR2,termTR3,termTR4)
             call cal_TermTR567 (EE,cU,cV,cW,termTR5,termTR6,termTR7)  
             call cal_TermTR9 (AA,cU,cV,cW,termTR9)
             str =  cU * nv * (nw+1) + (cV-1) * (nw+1) + cW + 1
             F(str) = termTR0_t + TermTR2-TermTR4-TermTR9-(TermTR5+TermTR6+TermTR7)/Lewis-TermTR1+TermTR3
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    s_min=ns ; s_max=1;
    !$OMP DO SCHEDULE(STATIC,(nu+1)/nb_taches)
    do cS=1,ns
       s_min=min(s_min,CS) ; s_max=max(s_max,CS)
       do cP=0,np
          do cT=0,nt
             call cal_TermEn0_t (dGGdt,cS,cP,cT,TermEn0_t)
             call cal_TermEn1234 (AA,BB,GG,cS,cP,cT,termE1,termE2,termE3,termE4)
             call cal_TermEn567 (GG,cS,cP,cT,termE5,termE6,termE7)  
             call cal_TermEn9 (BB,cS,cP,cT,termE9)
             sEn = nve + (cS-1) * (np+1) * (nt+1) + cP * (nt+1) + cT + 1
             F(sEn) =  TermEn0_t + TermE2 -TermE4  -TermE1+TermE3 -(TermE5+TermE6+TermE7) +TermE9
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    !$OMP END PARALLEL
    !$OMP END PARALLEL
    !$OMP END PARALLEL 

    !do i=1,nve+nvg
    !  norm = norm + F(i)**2.d0
    !  print*,F(i)
    !end do
    !norm = Dsqrt (norm)
    !print*,compt, norm

    return 
  end Subroutine Cal_res
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFX123 (cI,cJ,cK,termFx1,termFx2,termFx3)
    use in_data  
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK
    !real*8    , dimension (:,:)  , intent (in)  :: AA(0:ni,1:nj,1:nk)
    real*8                       , intent (out) :: termFx1,termFx2,termFx3
    real*8                                      :: cIr,cJr,cKr
    cIr = dble(cI)
    cJr = dble(cJ)
    cKr = dble(cK)

    termFX1 = -(pi**2.d0) *cIr**(2.d0)

    if (cI.eq.0) then  
       termFX2 = -2.d0* (pi**2.d0)  *cJr**(2.d0)
       termFX3 = -2.d0* (pi**2.d0)  *cKr**(2.d0)
    else 
       termFX2 = -(pi**2.d0) *cJr**(2.d0)
       termFX3 = -(pi**2.d0) *cKr**(2.d0)
    end if

  end subroutine cal_TermFX123
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFX4 (EE,cI,cJ,cK,termFx4)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK
    real*8    , dimension (:,:)  , intent (in)  :: EE(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: termFx4
    real*8                                      :: vr
    INTEGER :: v,w

    termFx4 = 0.d0 

    if (cI .le. nu) then 
       do v=1,nv
          do w=0,nw
             vr = dble(v)
             termFx4 = termFx4 + vr * matphi(cJ,v) * matphi(cK,w) * EE(cI,v,w)
          end do
       end do

       if (cI.eq.0) then 
          termFx4 = 2.d0*termFx4 * Ra / pi
       else 
          termFx4 = termFx4 * Ra /  pi
       end if
    end if
  end subroutine cal_TermFX4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFX5 (cI,cJ,cK,termFx5)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK
    real*8                       , intent (out) :: termFx5
    termFx5 = 0.d0 
    if (cI.eq.0)then
       termFx5 = 2*Ra*matphi(cJ,0)*matphi(cK,0)/ (pi**2.d0)
    else
       termFx5 = 0
    end if
  end subroutine cal_TermFX5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFX6 (GG,cI,cJ,cK,termFx6)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cI,cJ,cK
    real*8    , dimension (:,:)  , intent (in)  :: GG(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: termFx6
    real*8                                      :: pr
    INTEGER :: s,t

    termFx6 = 0.d0 

    if (cJ .le. np) then

       do s=1,ns
          do t=0,nt
             termFx6 = termFx6 - dble(cJ) * matphi(s,cI) * matphi(cK,t) * GG(s,cJ,t)
          end do

       end do
    end if


    termFx6 = termFx6 * Ra /  pi
  end subroutine cal_TermFX6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFY123 (cL,cM,cN,termFy1,termFy2,termFy3)
    use in_data  
    implicit none 
    integer                      , intent (in)  :: cL,cM,cN
    !real*8    , dimension (:,:)  , intent (in)  :: BB(1:nl,0:nm,1:nn)
    real*8                       , intent (out) :: termFy1,termFy2,termFy3
    real*8                                      :: cLr,cMr,cNr


    cLr = dble(cL)
    cMr = dble(cM)
    cNr = dble(cN)

    termFy2 = - (pi*cMr)**2.d0

    if (cM.eq.0) then  
       termFy1 = -2.d0* (pi**2.d0) *   cLr**(2.d0)
       termFy3 = -2.d0* (pi**2.d0) *   cNr**(2.d0)
    else 
       termFy1 = - (pi**2.d0) *  cLr**(2.d0)
       termFy3 = - (pi**2.d0) *  cNr**(2.d0)
    end if


  end subroutine cal_TermFY123
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFy4 (EE,cL,cM,cN,termFy4)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cL,cM,cN
    real*8    , dimension (:,:)  , intent (in)  :: EE(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: termFy4
    real*8                                      :: cLr
    INTEGER :: v,w

    termFy4 = 0.d0 
    cLr      = dble (cL)

    if (cL .le. nu) then 
       do v=1,nv
          do w=0,nw
             termFy4 = termFy4 +  matphi(v,cM) * matphi(cN,w) * EE(cL,v,w) 
          end do
       end do
    else 
       termFy4 = 0.d0 
    end if
    termFy4 = -cLr * Ra * termFy4  / pi


  end subroutine cal_TermFy4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFy5 (GG,cL,cM,cN,termFy5)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cL,cM,cN
    real*8    , dimension (:,:)  , intent (in)  :: GG(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: termFy5
    real*8                                      :: sr
    INTEGER :: s,t

    termFy5 = 0.d0 

    if (cM .le. np) then 
       do s=1,ns
          do t=0,nt
             sr = dble(s)
             termFy5 = termFy5 + sr * matphi(cL,s) * matphi(cN,t) * GG(s,cM,t)
          end do
       end do

       if (cM.eq.0) then 
          termFy5 = 2.d0*termFy5 * Ra / pi
       else 
          termFy5 = termFy5 * Ra /  pi
       end if
    end if
  end subroutine cal_TermFy5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermFy6 (cL,cM,cN,termFy6)
    use in_data  
    use matrices
    implicit none 
    integer                      , intent (in)  :: cL,cM,cN
    real*8                       , intent (out) :: termFy6
    termFy6 = 0.d0 
    if (cM.eq.0)then
       termFy6 = 2*Ra*matphi(cL,0)*matphi(cN,0)/ (pi**2.d0)
    else
       termFy6 = 0
    end if
  end subroutine cal_TermFy6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermTR0_t (dEEdt,cU,cV,cW,termTR0_t)
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8    , dimension (:,:)  , intent (in)  :: dEEdt(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: termTR0_t


    termTR0_t = poro * dEEdt (cU,cV,cW)

    if (Cu.eq. 0 .and. CW .eq. 0) termTR0_t = 4.d0 * poro * dEEdt (cU,cV,cW)
    if (Cu.eq. 0 .and. CW .ne. 0) termTR0_t = 2.d0 * poro * dEEdt (cU,cV,cW)
    if (Cu.ne. 0 .and. CW .eq. 0) termTR0_t = 2.d0 * poro * dEEdt (cU,cV,cW)

  end subroutine cal_TermTR0_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  cal_TermEn0_t (dGGdt,cS,cP,cT,TermEn0_t)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT
    real*8    , dimension (:,:)  , intent (in)  :: dGGdt(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: TermEn0_t


    TermEn0_t = sigma * dGGdt (cS,cP,cT)

    if (Cp.eq. 0 .and. CT .eq. 0) TermEn0_t = 4.d0 *  sigma * dGGdt (cS,cP,cT)
    if (Cp.eq. 0 .and. CT .ne. 0) TermEn0_t = 2.d0 *  sigma * dGGdt (cS,cP,cT)
    if (Cp.ne. 0 .and. CT .eq. 0) TermEn0_t = 2.d0 *  sigma * dGGdt (cS,cP,cT)

  end subroutine cal_TermEn0_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermTR1234 (AA,BB,EE,cU,cV,cW,termTR1,termTR2,termTR3,termTR4)
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8    , dimension (:,:)  , intent (in)  :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn),EE(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: termTR1,termTR2,termTR3,termTR4
    INTEGER :: i, j, k, l, m, n, u, v, w, ju, jv, jw

    real*8  :: rk, rv, rw, rj, rn, ru,rl
    termTR1=0.d0
    termTR2=0.d0
    termTR3=0.d0
    termTR4=0.d0

    do u=0,nu
       do v=1,nv
          do w=0,nw 

             rv = dble(v)
             rw = dble(w)
             ru = dble(u)


             if (cU.eq.0) then


                if ( cu-u .ge. 0 .and. cu-u .le. ni) then 
                   if (cv-v .ge. 1 .and. cv-v .le. nj) then
                      rj= dble(cV-v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,cv-v,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,cv-v,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,cv-v,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (cV+v .ge. 1 .and. cV+v .le. nj) then
                      rj= dble(cV+v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,cv+v,cW-w)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,cv+v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,cv+v,w-cW)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,cv+v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,cv+v,w+cw)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,cv+v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 1 .and. v-cV .le. nj) then 
                      rj= dble(v-cV)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2-(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,v-cV,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2-(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,v-cV,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2-(pi**2.d0)/4.0d0*rk*rv*AA(cu-u,v-cV,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(cu-u,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
                if ( u-cU .ge. 0 .and. u-cU .le. ni) then
                   if (cv-v .ge. 1 .and. cv-v .le. nj) then
                      rj= dble(cV-v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,cv-v,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,cv-v,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,cv-v,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (cV+v .ge. 1 .and. cV+v .le. nj) then
                      rj= dble(cV+v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,cv+v,cW-w)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,cv+v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,cv+v,w-cW)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,cv+v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,cv+v,w+cw)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,cv+v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 1 .and. v-cV .le. nj) then 
                      rj= dble(v-cV)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2-(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,v-cV,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2-(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,v-cV,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2-(pi**2.d0)/4.0d0*rk*rv*AA(u-cU,v-cV,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/4.0d0*rj*rw*AA(u-cU,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
             else
                if ( cu-u .ge. 0 .and. cu-u .le. ni) then 
                   if (cv-v .ge. 1 .and. cv-v .le. nj) then
                      rj= dble(cV-v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,cv-v,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,cv-v,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,cv-v,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (cV+v .ge. 1 .and. cV+v .le. nj) then
                      rj= dble(cV+v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,cv+v,cW-w)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,cv+v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,cv+v,w-cW)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,cv+v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,cv+v,w+cw)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,cv+v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 1 .and. v-cV .le. nj) then 
                      rj= dble(v-cV)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,v-cV,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,v-cV,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(cu-u,v-cV,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cu-u,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of d1
                if ( cU+u .ge. 0 .and. cU+u .le. ni) then 
                   if (cv-v .ge. 1 .and. cv-v .le. nj) then
                      rj= dble(cV-v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,cv-v,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,cv-v,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,cv-v,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (cV+v .ge. 1 .and. cV+v .le. nj) then
                      rj= dble(cV+v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,cv+v,cW-w)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,cv+v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,cv+v,w-cW)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,cv+v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,cv+v,w+cw)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,cv+v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 1 .and. v-cV .le. nj) then 
                      rj= dble(v-cV)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,v-cV,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,v-cV,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge.1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(cU+u ,v-cV,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(cU+u ,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of d2
                if ( u-cU .ge. 0 .and. u-cU .le. ni) then
                   if (cv-v .ge. 1 .and. cv-v .le. nj) then
                      rj= dble(cV-v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,cv-v,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,cv-v,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,cv-v,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (cV+v .ge. 1 .and. cV+v .le. nj) then
                      rj= dble(cV+v)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,cv+v,cW-w)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,cv+v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,cv+v,w-cW)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,cv+v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2+(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,cv+v,w+cw)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,cv+v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 1 .and. v-cV .le. nj) then 
                      rj= dble(v-cV)
                      if (cW-w .ge.1 .and. cW-w .le. nk) then              
                         rk = dble(cW-w)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,v-cV,cW-w)*EE(u,v,w)
                         termTR4=termTR4+(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nk) then
                         rk = dble(w-CW)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,v-cV,w-cW)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nk) then
                         rk = dble(w+CW)
                         termTR2=termTR2-(pi**2.d0)/8.0d0*rk*rv*AA(u-cU,v-cV,w+cw)*EE(u,v,w)
                         termTR4=termTR4-(pi**2.d0)/8.0d0*rj*rw*AA(u-cU,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if


!!!!!!!!!!!!!!!!!!!!!!!!!!end of d3

             end if

             !print*,termTR2,termTR4
             !pause 



             if (cU .eq. 0) then          !!!!!!!!!!!!!Starting TermTR 1&3
                if ( u-cU .ge. 1 .and. u-cU .le. nl) then
                   rl= dble(u-cU)
                   if (cV-v .ge. 0 .and. cV-v .le. nm) then
                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,cv-v,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rw*rl*BB(u-cU,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,cv-v,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(u-cU,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,cv-v,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(u-cU,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 0 .and. v-cV .le. nm) then

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,v-cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rw*rl*BB(u-cU,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,v-cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(u-cU,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,v-cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(u-cU,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v+cV .ge. 0 .and. v+cV .le. nm) then 

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,v+cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rl*rw*BB(u-cU,v+cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,v+cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rl*rw*BB(u-cU,v+cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(u-cU,v+cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rl*rw*BB(u-cU,v+cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
!!!!!!!!!!end of d14
                if ( cU-u .ge. 1 .and. cU-u .le. nl) then
                   rl= dble(cU-u)
                   if (cV-v .ge. 0 .and. cV-v .le. nm) then
                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,cv-v,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rw*rl*BB(cU-u,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,cv-v,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(cU-u,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,cv-v,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(cU-u,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 0 .and. v-cV .le. nm) then

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,v-cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(cU-u,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,v-cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(cU-u,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,v-cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rw*rl*BB(cU-u,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v+cV .ge.0 .and. v+cV .le. nm) then 

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,v+cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/4.0d0*rl*rw*BB(cU-u,v+cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,v+cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rl*rw*BB(cU-u,v+cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/4.0d0*rn*ru*BB(cU-u,v+cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/4.0d0*rl*rw*BB(cU-u,v+cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of d13
             else
                if ( cU-u .ge. 1 .and. cU-u .le. nl) then
                   rl= dble(cU-u)
                   if (cV-v .ge. 0 .and. cV-v .le. nm) then
                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,cv-v,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rw*rl*BB(cU-u,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,cv-v,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(cU-u,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,cv-v,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(cU-u,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 0 .and. v-cV .le. nm) then

                      if (cW-w .ge. 1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,v-cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rw*rl*BB(cU-u,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,v-cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(cU-u,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,v-cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(cU-u,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v+cV .ge. 0 .and. v+cV .le. nm) then 

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,v+cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rl*rw*BB(cU-u,v+cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,v+cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rl*rw*BB(cU-u,v+cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(cU-u,v+cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rl*rw*BB(cU-u,v+cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
                if ( u-cU .ge. 1 .and. u-cU .le. nl) then
                   rl= dble(u-cU)
                   if (cV-v .ge. 0 .and. cV-v .le. nm) then
                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,cv-v,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rw*rl*BB(u-cU,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,cv-v,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u-cU,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,cv-v,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u-cU,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge.0 .and. v-cV .le. nm) then

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,v-cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rw*rl*BB(u-cU,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,v-cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u-cU,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge.1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,v-cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u-cU,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v+cV .ge. 0 .and. v+cV .le. nm) then 

                      if (cW-w .ge.1.and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,v+cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rl*rw*BB(u-cU,v+cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,v+cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rl*rw*BB(u-cU,v+cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(u-cU,v+cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rl*rw*BB(u-cU,v+cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
                if ( u+cU .ge. 1 .and. u+cU .le. nl) then
                   rl= dble(u+cU)
                   if (cV-v .ge. 0 .and. cV-v .le. nm) then
                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,cv-v,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rw*rl*BB(u+cU,cv-v,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,cv-v,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u+cU,cv-v,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,cv-v,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u+cU,cv-v,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v-cV .ge. 0 .and. v-cV .le. nm) then

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,v-cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rw*rl*BB(u+cU,v-cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,v-cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u+cU,v-cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1-(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,v-cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rw*rl*BB(u+cU,v-cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                   if (v+cV .ge. 0 .and. v+cV .le. nm) then 

                      if (cW-w .ge.1 .and. cW-w .le. nn) then              
                         rn = dble(cW-w)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,v+cV,cW-w)*EE(u,v,w)
                         termTR3=termTR3-(pi**2.d0)/8.0d0*rl*rw*BB(u+cU,v+cV,cW-w)*EE(u,v,w)
                      end if
                      if ( w-cW .ge. 1 .and. w-cW .le. nn) then
                         rn = dble(w-CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,v+cV,w-cW)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rl*rw*BB(u+cU,v+cV,w-cW)*EE(u,v,w)
                      end if
                      if ( w+CW .ge. 1 .and. w+CW .le. nn) then
                         rn = dble(w+CW)
                         termTR1=termTR1+(pi**2.d0)/8.0d0*rn*ru*BB(u+cU,v+cV,w+cw)*EE(u,v,w)
                         termTR3=termTR3+(pi**2.d0)/8.0d0*rl*rw*BB(u+cU,v+cV,w+cW)*EE(u,v,w)
                      end if
                   end if
                end if
             end if
!!!!!!

          end do
       end do
    end do



  end subroutine cal_TermTR1234
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermTR567 (EE,cU,cV,cW,termTR5,termTR6,termTR7)       
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8    , dimension (:,:)  , intent (in)  :: EE(0:nu,1:nv,0:nw)
    real*8                       , intent (out) :: termTR5,termTR6,termTR7
    real*8                                      :: cUr,cVr,cWr

    cUr=dble(cU)
    cVr=dble(cV)
    cWr=dble(cW)

    if (cW.eq.0) then
       termTR5=-2.d0*(pi**2.d0)*(cUr**2.d0)*EE(cU,cV,cW)
    else
       termTR5=-(pi**2.d0)*(cUr**2.d0)*EE(cU,cV,cW) 
    endif

    if ((cW.eq.0).and.(cU.eq.0)) then
       termTR6=-4.d0*(pi**2.d0)*(cVr**2.d0)*EE(cU,cV,cW) 
    elseif ((cW.ne.0).and.(cU.ne.0)) then
       termTR6=-(pi**2.d0)*(cVr**2.d0)*EE(cU,cV,cW)  
    else
       termTR6=-2.d0*(pi**2.d0)*(cVr**2.d0)*EE(cU,cV,cW)
    endif

    if (cU.eq.0) then
       termTR7=-2.d0*(pi**2.d0)*(cWr**2.d0)*EE(cU,cV,cW)
    else
       termTR7=-(pi**2.d0)*(cWr**2.d0)*EE(cU,cV,cW) 
    endif

  end subroutine cal_TermTR567
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermTR8 (BB,cU,cV,cW,termTR8)       
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8    , dimension (:,:)  , intent (in)  :: BB(1:nl,0:nm,1:nn)
    real*8                       , intent (out) :: termTR8
    real*8                                      :: t1,t2
    INTEGER :: l,m

    termTR8=0.d0
    if ((cW.le.nn).and.(cW.ne.0)) then
       if (cU.eq.0) then
          do m=0,nm
             t1=dble(cV+m)
             t2=dble(cV-m)
             if (m.ne.cV) then
                termTR8=termTR8+pi/2.d0*dble(cW)*(1.D0/t1+1.D0/t2)*BB(1,m,cW)
             else
                termTR8=termTR8+pi/4.d0*dble(cW)/dble(cV)*BB(1,cV,cW)
             endif
          enddo
       else 
          do l=1,nl
             if ((l.eq.(1-cU)).or.(l.eq.(1+cU))) then
                do m=0,nm
                   t1=dble(cV+m)
                   t2=dble(cV-m)
                   if (m.ne.cV) then
                      termTR8=termTR8+pi/4.d0*dble(cW)*(1.D0/t1+1.D0/t2)*BB(l,m,cW)
                   else
                      termTR8=termTR8+pi/8.d0*dble(cW)/dble(cV)*BB(l,cV,cW)
                   endif
                enddo
             elseif (l.eq.(cU-1)) then   
                do m=0,nm
                   t1=dble(cV+m)
                   t2=dble(cV-m)
                   if (m.ne.cV) then
                      termTR8=termTR8-pi/4.d0*dble(cW)*(1.D0/t1+1.D0/t2)*BB(l,m,cW)
                   else
                      termTR8=termTR8-pi/8.d0*dble(cW)/dble(cV)*BB(l,cV,cW)
                   endif
                enddo
             endif
          enddo
       endif
    endif

  end subroutine cal_TermTR8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermTR9 (AA,cU,cV,cW,termTR9)       
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8    , dimension (:,:)  , intent (in)  :: AA(0:ni,1:nj,1:nk)
    real*8                       , intent (out) :: termTR9
    INTEGER :: i,j,k

    termTR9=0.0
    if (cW.ne.0) then
       if (cU.eq.0) then
          termTR9=pi*dble(cW)*AA(0,cV,cW)    !When i can be not equal to u but u=0 can add some values to the Term9
       end if
    end if

    do i=0,ni
       do j=1,nj
          do k=1,nk
             if ((i.eq.cU).and. (j.eq.cV).and.(k.eq.cW)) then
                termTR9=termTR9+pi*dble(k)*AA(i,j,k)
             end if
          enddo
       enddo
    enddo

  end subroutine cal_TermTR9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermTR10 (cU,cV,cW,termTR10)       
    use in_data
    implicit none
    integer                      , intent (in)  :: cU,cV,cW
    real*8                       , intent (out) :: termTR10

    termTR10=0.d0
    if ((cU.eq.1).and.(cW.eq.0)) termTR10=2.D0*pi/dble(cV)

  end subroutine cal_TermTR10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermEn1234(AA,BB,GG,cS,cP,cT,termEn1,termEn2,termEn3,termEn4)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT
    real*8    , dimension (:,:)  , intent (in)  :: AA(0:ni,1:nj,1:nk),BB(1:nl,0:nm,1:nn),GG(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: termEn1,termEn2,termEn3,termEn4
    INTEGER :: i, j, k, l, m, n, s, p, t, ju, jv, jw

    real*8  :: rp, rt,rs, rj, rk,rn,rl
    termEn1=0.d0
    termEn2=0.d0
    termEn3=0.d0
    termEn4=0.d0
    do s=1,ns
       do p=0,np
          do t=0,nt
             rp=dble(p)
             rt=dble(t)
             rs=dble(s)
             if ( cS-s.ge.0 .and. cS-s.le.ni) then
                if (cP+p.ge.1 .and. cP+p.le.nj) then
                   rj=dble(cP+p)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,cP+p,cT-t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,cP+p,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,cP+p,cT+t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,cP+p,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,cP+p,t-cT)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,cP+p,t-cT)*GG(s,p,t)
                   end if
                end if
                if (p-cP.ge.1 .and. p-cP.le.nj) then
                   rj=dble(p-cP)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,p-cP,cT-t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,p-cP,cT-t)*GG(s,p,t)

                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,p-cP,cT+t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,p-cP,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,p-cP,t-cT)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,p-cP,t-cT)*GG(s,p,t)
                   end if
                end if
                if (cP-p.ge.1 .and. cP-p.le.nj) then
                   rj=dble(cP-p)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,cP-p,cT-t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,cP-p,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,cP-p,cT+t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,cP-p,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(cS-s,cP-p,t-cT)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(cS-s,cP-p,t-cT)*GG(s,p,t)
                   end if
                end if
             end if
             if ( s-cS.ge.0 .and. s-cS.le.ni) then
                if (cP+p.ge.1 .and. cP+p.le.nj) then
                   rj=dble(cP+p)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,cP+p,cT-t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,cP+p,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,cP+p,cT+t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,cP+p,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,cP+p,t-cT)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,cP+p,t-cT)*GG(s,p,t)
                   end if
                end if
                if (p-cP.ge.1 .and. p-cP.le.nj) then
                   rj=dble(p-cP)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,p-cP,cT-t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,p-cP,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,p-cP,cT+t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,p-cP,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,p-cP,t-cT)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,p-cP,t-cT)*GG(s,p,t)
                   end if
                end if
                if (cP-p.ge.1 .and. cP-p.le.nj) then
                   rj=dble(cP-p)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,cP-p,cT-t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,cP-p,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,cP-p,cT+t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,cP-p,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s-cS,cP-p,t-cT)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s-cS,cP-p,t-cT)*GG(s,p,t)
                   end if
                end if
             end if
             if ( s+cS.ge.0 .and. s+cS.le.ni) then
                if (cP+p.ge.1 .and. cP+p.le.nj) then
                   rj=dble(cP+p)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,cP+p,cT-t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,cP+p,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,cP+p,cT+t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,cP+p,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,cP+p,t-cT)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,cP+p,t-cT)*GG(s,p,t)
                   end if
                end if
                if (p-cP.ge.1 .and. p-cP.le.nj) then
                   rj=dble(p-cP)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,p-cP,cT-t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,p-cP,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,p-cP,cT+t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,p-cP,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2+(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,p-cP,t-cT)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,p-cP,t-cT)*GG(s,p,t)
                   end if
                end if
                if (cP-p.ge.1 .and. cP-p.le.nj) then
                   rj=dble(cP-p)
                   if (cT-t .ge. 1 .and. cT-t .le. nk) then
                      rk=dble(cT-t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,cP-p,cT-t)*GG(s,p,t)
                      termEn4=termEn4-(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,cP-p,cT-t)*GG(s,p,t)
                   end if
                   if (cT+t .ge. 1 .and. cT+t .le. nk) then
                      rk=dble(cT+t)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,cP-p,cT+t)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,cP-p,cT+t)*GG(s,p,t)
                   end if
                   if (t-cT .ge. 1 .and. t-cT .le. nk) then
                      rk=dble(t-cT)
                      termEn2=termEn2-(pi**2.d0)/8.0d0*rk*rp*AA(s+cS,cP-p,t-cT)*GG(s,p,t)
                      termEn4=termEn4+(pi**2.d0)/8.0d0*rj*rt*AA(s+cS,cP-p,t-cT)*GG(s,p,t)
                   end if
                end if
             end if
             !starting term 1,3 
             if (cP .eq. 0) then
                if (cS-s .ge. 1 .and. cS-s .le. nl) then
                   rl=dble(cS-s)
                   if (cP-p .ge. 0 .and. cP-p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS-s,cP-p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(cS-s,cP-p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS-s,cP-p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(cS-s,cP-p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS-s,cP-p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(cS-s,cP-p,t-cT)*GG(s,p,t) 
                      end if
                   end if

                   if (cP+p .ge. 0 .and. cP+p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS-s,cP+p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(cS-s,cP+p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS-s,cP+p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(cS-s,cP+p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS-s,cP+p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(cS-s,cP+p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                end if
                if (cS+s .ge. 1 .and. cS+s .le. nl) then
                   rl=dble(cS+s)
                   if (cP-p .ge. 0 .and. cP-p  .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS+s,cP-p ,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(cS+s,cP-p ,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS+s,cP-p ,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(cS+s,cP-p ,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS+s,cP-p ,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(cS+s,cP-p ,t-cT)*GG(s,p,t) 
                      end if
                   end if

                   if (p+cP .ge. 0 .and. p+cP.le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS+s,p+cP,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(cS+s,p+cP,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS+s,p+cP,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(cS+s,p+cP,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/4.0d0*rn*rs*BB(cS+s,p+cP,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(cS+s,p+cP,t-cT)*GG(s,p,t) 
                      end if
                   end if
                end if
                if (s-cS .ge. 1 .and. s-cS .le. nl) then
                   rl=dble(s-cS)
                   if (cP-p .ge. 0 .and. cP-p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1-(pi**2.d0)/4.0d0*rn*rs*BB(s-cS,cP-p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(s-cS,cP-p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1-(pi**2.d0)/4.0d0*rn*rs*BB(s-cS,cP-p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(s-cS,cP-p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1-(pi**2.d0)/4.0d0*rn*rs*BB(s-cS,cP-p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(s-cS,cP-p,t-cT)*GG(s,p,t) 
                      end if
                   end if

                   if (cP+p .ge. 0 .and. cP+p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1-(pi**2.d0)/4.0d0*rn*rs*BB(s-cS,cP+p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/4.0d0*rl*rt*BB(s-cS,cP+p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1-(pi**2.d0)/4.0d0*rn*rs*BB(s-cS,cP+p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(s-cS,cP+p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1-(pi**2.d0)/4.0d0*rn*rs*BB(s-cS,cP+p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/4.0d0*rl*rt*BB(s-cS,cP+p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                end if
             else
                if (cS-s .ge. 1 .and. cS-s .le. nl) then
                   rl=dble(cS-s)
                   if (cP-p .ge. 0 .and. cP-p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,cP-p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,cP-p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,cP-p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,cP-p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,cP-p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,cP-p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                   if (cP+p .ge. 0 .and. cP+p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,cP+p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,cP+p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,cP+p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,cP+p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,cP+p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,cP+p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                   if (p-cP .ge. 0 .and. p-cP .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,p-cP,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,p-cP,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,p-cP,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,p-cP,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS-s,p-cP,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS-s,p-cP,t-cT)*GG(s,p,t) 
                      end if
                   end if
                end if
                if (cS+s .ge. 1 .and. cS+s .le. nl) then
                   rl=dble(cS+s)
                   if (cP-p .ge. 0 .and. cP-p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,cP-p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,cP-p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,cP-p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,cP-p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,cP-p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,cP-p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                   if (cP+p .ge. 0 .and. cP+p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,cP+p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,cP+p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,cP+p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,cP+p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,cP+p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,cP+p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                   if (p-cP .ge. 0 .and. p-cP .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,p-cP,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,p-cP,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,p-cP,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,p-cP,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1+(pi**2.d0)/8.0d0*rn*rs*BB(cS+s,p-cP,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(cS+s,p-cP,t-cT)*GG(s,p,t) 
                      end if
                   end if
                end if
                if (s-cS .ge. 1 .and. s-cS .le. nl) then
                   rl=dble(s-cS)
                   if (cP-p .ge. 0 .and. cP-p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,cP-p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,cP-p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,cP-p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,cP-p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,cP-p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,cP-p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                   if (cP+p .ge. 0 .and. cP+p .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,cP+p,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,cP+p,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,cP+p,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,cP+p,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,cP+p,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,cP+p,t-cT)*GG(s,p,t) 
                      end if
                   end if
                   if (p-cP .ge. 0 .and. p-cP .le. nm) then
                      if ( cT-t .ge. 1 .and. cT-t .le.nn) then
                         rn=dble(cT-t)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,p-cP,cT-t)*GG(s,p,t)
                         TermEn3=TermEn3+(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,p-cP,cT-t)*GG(s,p,t)
                      end if
                      if (cT+t .ge. 1 .and. cT+t .le. nn) then
                         rn=dble(cT+t)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,p-cP,cT+t)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,p-cP,cT+t)*GG(s,p,t) 
                      end if
                      if (t-cT .ge. 1 .and. t-cT .le. nn) then
                         rn=dble(t-cT)
                         TermEn1=TermEn1-(pi**2.d0)/8.0d0*rn*rs*BB(s-cS,p-cP,t-cT)*GG(s,p,t)
                         TermEn3=TermEn3-(pi**2.d0)/8.0d0*rl*rt*BB(s-cS,p-cP,t-cT)*GG(s,p,t) 
                      end if
                   end if
                end if
             end if

          end do
       end do
    end do
  end subroutine cal_TermEn1234

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermEn567(GG,cS,cP,cT,termE5,termE6,termE7)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT
    real*8    , dimension (:,:)  , intent (in)  :: GG(1:ns,0:np,0:nt)
    real*8                       , intent (out) :: termE5,termE6,termE7
    real*8                                      :: cSr,cPr,cTr

    cSr=dble(cS)
    cPr=dble(cP)
    cTr=dble(cT)

    if ((cT .eq. 0) .and. (cP .eq. 0)) then
       TermE5=-4.d0*(pi**2.d0)*(cSr**2.d0)*GG(cS,cP,cT)
    elseif ((cT.ne.0).and.(cP.ne.0)) then
       termE5=-(pi**2.d0)*(cSr**2.d0)*GG(cS,cP,cT)  
    else
       termE5=-2.d0*(pi**2.d0)*(cSr**2.d0)*GG(cS,cP,cT)
    endif

    if (cT.eq.0) then
       termE6=-2.d0*(pi**2.d0)*(cPr**2.d0)*GG(cS,cP,cT)
    else
       termE6=-(pi**2.d0)*(cPr**2.d0)*GG(cS,cP,cT) 
    endif
    !

    !
    if (cP.eq.0) then
       termE7=-2.d0*(pi**2.d0)*(cTr**2.d0)*GG(cS,cP,cT)
    else
       termE7=-(pi**2.d0)*(cTr**2.d0)*GG(cS,cP,cT) 
    endif

  end subroutine cal_TermEn567
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cal_TermEn9 (BB,cS,cP,cT,termE9)
    use in_data
    implicit none
    integer                      , intent (in)  :: cS,cP,cT
    real*8    , dimension (:,:)  , intent (in)  :: BB(1:nl,0:nm,1:nn)
    real*8                       , intent (out) :: termE9
    INTEGER :: l,m,n

    termE9=0.0
    if (cT.ne.0) then
       if (cP.eq.0) then
          termE9=pi*dble(cT)*BB(cS,0,cT)    
       end if
    end if
    do l=1,nl
       do m=0,nm
          do n=1,nn


             if ((l.eq.cS).and. (m.eq.cP).and.(n.eq.cT)) then

                termE9=termE9+pi*dble(n)*BB(l,m,n)

             end if

          enddo
       enddo
    enddo

  end subroutine cal_TermEn9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Calc_AA(EE,GG,AA)
    use in_data
    implicit none 
    real*8, dimension(:,:), intent(in) :: EE(0:nu, 1:nv, 0:nw), GG(1:ns,0:np,0:nt)
    integer :: CI, CJ, CK
    real*8, dimension(:,:), intent(out) :: AA(0:ni, 1:nj, 1:nk)
    real*8                              :: TermFx1,termFx2,termFx3,TermFx4,termFx5,termFx6
    integer                             :: i,j,k

    do CI=0,ni
       do CJ=1,nj
          do Ck=1,nk
             call cal_TermFX4 (EE,cI,cJ,cK,termFx4)
             call cal_TermFX5 (cI,cJ,cK,termFx5)
             call cal_TermFX6 (GG,cI,cJ,cK,termFx6)
             call cal_TermFX123 (cI,cJ,Ck,TermFx1,termFx2,termFx3)
             AA(ci,cj,ck)=(Ngrav*TermFx4 -Ngrav*TermFx5 +TermFx6 )/ (TermFx1+TermFx2+TermFx3)
          end do
       end do
    end do
  end subroutine Calc_AA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Calc_BB(EE,GG,BB)
    use in_data
    implicit none 
    real*8, dimension(:,:), intent(in) :: EE(0:nu, 1:nv, 0:nw), GG(1:ns,0:np,0:nt)
    integer:: CL, CM, CN
    real*8, dimension(:,:), intent(out) :: BB(1:nl, 0:nm, 1:nn)
    real*8                              :: TermFy1,termFy2,termFy3,TermFy4,termFy5,termFy6
    integer                             :: l,m,n

    do cl=1,nl
       do cm=0,nm
          do cn=1,nn

             call cal_TermFY123 (cL,cM,cN,TermFy1,termFy2,termFy3)
             call cal_TermFy4 (EE,cL,cM,cN,termFy4)
             call cal_TermFy5 (GG,cL,cM,cN,termFy5)
             call cal_TermFy6 (cL,cM,cN,termFy6)
             BB(cl,cm,cn)=(-Ngrav * TermFy4 - TermFy5 +TermFy6)/(TermFy1 + TermFy2 +TermFy3 )          
          end do
       end do
    end do
  end subroutine Calc_BB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end MODULE residual_module
